
/**
 * @file branch2.c
 *
 * @brief libsea3 API implementations
 *
 * @author Hajime Suzuki
 * @date 2016/1/11
 * @license Apache v2
 */
#include <stdint.h>
#include "sea.h"
#include "util/util.h"
#include "arch/arch.h"

/* aliasing vector macros */
#define _VECTOR_ALIAS_PREFIX	v32i8
#include "arch/vector_alias.h"

#define BW 				( 32 )
#define BLK 			( 32 )
#define MIN_BULK_BLOCKS	( 32 )
#define INIT_STACK_SIZE			( (uint64_t)32 * 1024 * 1024 )

/**
 * @macro _match
 * @brief alias to sequence matcher macro
 */
#ifndef _match
// #define _match			_eq 		/* for ASCII encoded sequence */
#define _match		_or 		/* for 2bit encoded */
// #define _match		_and		/* for 4bit encoded */
#endif /* _match */

/**
 * @macro _likely, _unlikely
 * @brief branch prediction hint for gcc-compatible compilers
 */
#define _likely(x)		__builtin_expect((x), 1)
#define _unlikely(x)	__builtin_expect((x), 0)

/**
 * @macro _force_inline
 * @brief inline directive for gcc-compatible compilers
 */
#define _force_inline	inline

/**
 * @macro _load_sc
 * @brief load constant vectors from sea_dp_context_s
 */
#define _load_sc(this, name)	( _bc_v16i8(_load_v16i8((this)->scv.name)) )

/* direction macros */
#define DYNAMIC 		( 1 )
#if DYNAMIC
#define direction_prefix 		dynamic_

/* direction determiners for the dynamic band algorithms */
/**
 * @macro _dir_update
 * @brief update direction determiner for the next band
 */
#define _dir_update(_dir, _vector) { \
	(_dir).dynamic.acc += (_ext(_vector, 0) - _ext(_vector, BW-1)); \
	(_dir).dynamic.array <<= 1; \
	(_dir).dynamic.array |= ((_dir).dynamic.acc < 0); \
}
/**
 * @macro _dir_adjust_reminder
 * @brief adjust direction array when termination is detected in the middle of the block
 */
#define _dir_adjust_reminder(_dir, i) { \
	(_dir).dynamic.array <<= (BLK - (i) - 1); \
}
/**
 * @macro _dir_test_bound, _dir_test_bound_cap
 * @brief test if the bound of the direction array is invaded
 */
#define _dir_test_bound(_dir, k, p) 		( /* nothing to do */ )
#define _dir_test_bound_cap(_dir, k, p)		( /* nothing to do */ )
/**
 * @macro _dir_is_down, _dir_is_right
 * @brief direction indicator (_dir_is_down returns true if dir == down)
 */
#define _dir_is_down(_dir)					( (_dir).dynamic.array & 0x01 )
#define _dir_is_right(_dir)					( !_dir_is_down(_dir) )

#else

#define direction_prefix 		guided_

/* direction determiners for the guided band algorithms */
#define _dir_update(_dir, _vector) { \
	(_dir).guided.ptr++; \
}
#define _dir_adjust_reminder(_dir, i)		{ /* nothing to do */ }
#define _dir_test_bound(_dir, k, p) ( \
	((k)->tdr - (_dir).guided.ptr + 2) - (p) - BLK \
)
#define _dir_test_bound_cap(_dir, k, p) ( \
	((k)->tdr - (_dir).guided.ptr + 2) - (p) \
)
#define _dir_is_down(_dir)					( *(_dir).guided.ptr != 0 )
#define _dir_is_right(_dir)					( *(_dir).guided.ptr == 0 )

#endif

/**
 * seqreader macros
 */
#define _rd_pos1(k)		_pv_v2i64(&(k)->rr.s.body.apos)
#define _rd_len1(k)		_pv_v2i64(&(k)->rr.s.body.alen)
#define _rd_pos2(k)		_pv_v2i64(&(k)->rr.s.tail.apos)
#define _rd_len2(k)		_pv_v2i64(&(k)->rr.s.tail.alen)
#define _rd_cnt(k)		_pv_v2i64(&(k)->rr.acnt)
// #define _rd_len(k)		_pv_v2i64(&(k)->work.lena)
#define _rd_bufa_base(k)		( (k)->rr.bufa + BLK + BW )
#define _rd_bufb_base(k)		( (k)->rr.bufb )
#define _rd_bufa(k, pos, len)	( _rd_bufa_base(k) - (pos) - (len) )
#define _rd_bufb(k, pos, len)	( _rd_bufb_base(k) + (pos) )
#define _lo64(v)		_ext_v2i64(v, 0)
#define _hi64(v)		_ext_v2i64(v, 1)

/**
 * @fn rd_load_section
 * @brief store given sea_section_pair_s to current working buffer on the dp_context
 */
static _force_inline
void rd_load_section(
	struct sea_dp_context_s *this,
	struct sea_section_pair_s *ptr_sec)
{
	_memcpy_blk_au(&this->rr.s, ptr_sec, sizeof(struct sea_section_pair_s));
	return;
}

/**
 * @fn rd_save_section
 * @brief save current section pair on the context to ptr_sec
 */
static _force_inline
void rd_save_section(
	struct sea_dp_context_s *this,
	struct sea_section_pair_s *ptr_sec)
{
	_memcpy_blk_ua(ptr_sec, &this->rr.s, sizeof(struct sea_section_pair_s));
	return;
}

/**
 * @fn rd_load_seq
 * @brief load previous fetched sequence onto seq buffer on the context
 */
static _force_inline
void rd_load_seq(
	struct sea_dp_context_s *this,
	struct sea_joint_tail_s *tail)
{
	_store(_rd_bufa(this, 0, BW), _load(tail->wa));
	_store(_rd_bufb(this, 0, BW), _load(tail->wb));
	return;
}

/**
 * @fn rd_save_seq
 * @brief save current seq buffer to joint_tail
 */
static _force_inline
void rd_save_seq(
	struct sea_dp_context_s *this,
	struct sea_joint_tail_s *tail)
{
	_store(tail->wa, _load(_rd_bufa(this, 0, BW)));
	_store(tail->wb, _load(_rd_bufb(this, 0, BW)));
	return;
}

/**
 * @fn rd_go_down, rd_go_right
 * @brief increment counter
 */
static _force_inline
void rd_go_right(
	struct sea_dp_context_s *this)
{
	this->rr.acnt++;
	return;
}
static _force_inline
void rd_go_down(
	struct sea_dp_context_s *this)
{
	this->rr.bcnt++;
	return;
}

/**
 * @fn rd_bulk_fetch
 * @brief fetch sequence to seq buffer (fast)
 */
static _force_inline
void rd_bulk_fetch(
	struct sea_dp_context_s *this)
{
	/* load pos, len and cnt */
	v2i64_t pos = _load_v2i64(_rd_pos1(this));
	v2i64_t len = _load_v2i64(_rd_len1(this));
	v2i64_t cnt = _load_v2i64(_rd_cnt(this));

	/* update pos and len */
	pos = _add_v2i64(pos, cnt);
	len = _add_v2i64(len, cnt);
	_store_v2i64(_rd_pos1(this), pos);
	_store_v2i64(_rd_len1(this), len);

	/* fetch seq a */
	vec_t t = _loadu(_rd_bufa(this, _lo64(cnt), BW));
	_rd_load(this->r.loada,
		_rd_bufa(this, BW, BLK),			/* dst */
		this->rr.p.pa,						/* src */
		rev(_lo64(pos), this->rr.p.alen),	/* pos */
		this->rr.p.alen,					/* lim */
		BLK);								/* len */
	_store(_rd_bufa(this, 0, BW), t);

	/* fetch seq b */
	_store(_rd_bufb(this, 0, BW), _loadu(_rd_bufb(this, _hi64(cnt), BW)));
	_rd_load(this->r.loadb,
		_rd_bufb(this, BW, BLK),			/* dst */
		this->rr.p.pb,						/* src */
		_hi64(pos),							/* pos */
		this->rr.p.blen,					/* lim */
		BLK);								/* len */

	/* clear counter */
	_store_v2i64(_rd_cnt(this), _zero_v2i64());

	return;
}

/**
 * @fn rd_test_bulk_fetch
 * @brief check if bulk fetch is available
 */
static _force_inline
int64_t rd_test_fast_fetch(
	struct sea_dp_context_s *this,
	uint32_t p)
{
	return(((int64_t)this->rr.s.body.alen - this->rr.acnt - BW)
		 | ((int64_t)this->rr.s.body.blen - this->rr.bcnt - BW)
		 | ((int64_t)this->rr.s.limp - p));
}

/**
 * @fn rd_cap_fetch
 * @brief fetch sequence to seq buffer (for cap fill)
 */

static _force_inline
void rd_cap_fetch(
	struct sea_dp_context_s *this)
{
	/* constants */
	v2i64_t const tot = _set_v2i64(BLK);
	v2i64_t const zero = _zero_v2i64();

	/* load lengths */
	v2i64_t len1 = _load_v2i64(_rd_len1(this));
	v2i64_t len2 = _load_v2i64(_rd_len2(this));

	/* load cnt */
	v2i64_t cnt = _load_v2i64(_rd_cnt(this));
	v2i64_t cnt2 = _max_v2i64(_sub_v2i64(cnt, len1), zero);
	v2i64_t cnt1 = _sub_v2i64(cnt, cnt2);

	/* update section 1 */
	v2i64_t pos1 = _add_v2i64(_load_v2i64(_rd_pos1(this)), cnt1);
	len1 = _min_v2i64(_sub_v2i64(len1, cnt1), tot);
	_store_v2i64(_rd_pos1(this), pos1);
	_store_v2i64(_rd_len1(this), len1);

	/* update section 2 */
	v2i64_t pos2 = _add_v2i64(_load_v2i64(_rd_pos2(this)), cnt2);
	len2 = _min_v2i64(_sub_v2i64(len2, cnt2), _sub_v2i64(tot, len1));
	_store_v2i64(_rd_pos2(this), pos2);
	_store_v2i64(_rd_len2(this), len2);

	/* fetch seq a */
	vec_t t = _loadu(_rd_bufa(this, _lo64(cnt), BW));
	_rd_load(this->r.loada,
		_rd_bufa(this, BW + _lo64(len1), _lo64(len2)),	/* dst */
		this->rr.p.pa,						/* src */
		rev(_lo64(pos2), this->rr.p.alen),	/* pos */
		this->rr.p.alen,					/* lim */
		_lo64(len2));						/* len */
	_rd_load(this->r.loada,
		_rd_bufa(this, BW, _lo64(len1)),	/* dst */
		this->rr.p.pa,						/* src */
		rev(_lo64(pos2), this->rr.p.alen),	/* pos */
		this->rr.p.alen,					/* lim */
		_lo64(len1));						/* len */
	_store(_rd_bufa(this, 0, BW), t);

	/* fetch seq b */
	_store(_rd_bufb(this, 0, BW), _loadu(_rd_bufb(this, _hi64(cnt), BW)));
	_rd_load(this->r.loadb,
		_rd_bufb(this, BW, _hi64(len1)),	/* dst */
		this->rr.p.pb,						/* src */
		_hi64(pos1),						/* pos */
		this->rr.p.blen,					/* lim */
		_hi64(len1));						/* len */
	_rd_load(this->r.loadb,
		_rd_bufb(this, BW + _hi64(len1), _hi64(len2)),
		this->rr.p.pb,						/* src */
		_hi64(pos2),						/* pos */
		this->rr.p.blen,					/* lim */
		_hi64(len2));						/* len */

	/* clear counter */
	_store_v2i64(_rd_cnt(this), zero);

	return;
}

/**
 * fill-in macros
 */
/**
 * @macro _fill_load_contexts
 * @brief load vectors onto registers
 */
#define _fill_load_contexts(_blk) \
	/* load direction determiner */ \
	union sea_dir_u dir = (_blk)->dir; \
	/* load large offset */ \
	int64_t offset = (_blk)->offset; \
	/* load vector registers */ \
	vec_t dh = _load(_pv((_blk)->diff.dh)); \
	vec_t dv = _load(_pv((_blk)->diff.dv)); \
	vec_t de = _load(_pv((_blk)->diff.de)); \
	vec_t df = _load(_pv((_blk)->diff.df)); \
	/* load delta vectors */ \
	vec_t delta = _load(_pv((_blk)->sd.delta)); \
	vec_t max = _load(_pv((_blk)->sd.max)); \
/**
 * @macro _fill_body
 * @brief update vectors
 */
#define _fill_body() { \
	vec_t _t = _match( \
		_loadu(_rd_bufa(this, 0, BW)), \
		_loadu(_rd_bufb(this, 0, BW))); \
	_t = _shuf(_t, _load_sc(this, sbv)); \
	_t = _max(_t, de); _t = _max(_t, df); \
	de = _max(de, dv); df = _max(df, dh); \
	vec_t _dh = _sub(_t, dv); \
	vec_t _dv = _sub(_t, dh); \
	vec_t _de = _sub(de, dh); \
	vec_t _df = _sub(df, dv); \
	*mask_ptr++ = (struct sea_mask_pair_s) { \
		_mask(_eq(_dh, _df)), _mask(_eq(_dv, _de)) \
	}; \
	dh = _dh; dv = _dv; \
	de = _add(_de, _load_sc(this, geav)); \
	df = _add(_df, _load_sc(this, gebv)); \
}
/**
 * @macro _fill_update_delta
 * @brief update small delta vector and max vector
 */
#define _fill_update_delta(_vector, _offset) { \
	delta = _add(delta, _vector); \
	delta = _add(delta, _offset); \
	max = _max(max, delta); \
	_dir_update(dir, _vector); \
}
/**
 * @macro _fill_right, _fill_down
 * @brief wrapper of _fill_body and _fill_update_delta
 */
#define _fill_right() { \
	rd_go_right(this);	/* increment sequence buffer pointer */ \
	dh = _shl(dh, 1);	/* shift left dh */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(dh, _load_sc(this, giav)); \
}
#define _fill_down() { \
	rd_go_down(this);	/* increment sequence buffer pointer */ \
	dv = _shr(dv, 1);	/* shift right dv */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(dv, _load_sc(this, gibv)); \
}
/**
 * @macro _fill_update_offset
 * @brief update offset and max vector, reset the small delta
 */
#define _fill_update_offset() { \
	int8_t _cd = _ext(delta, BW/2); \
	offset += _cd; \
	delta = _sub(delta, _set(_cd)); \
	max = _sub(max, _set(_cd)); \
}
/**
 * @macro _fill_store_vectors
 * @brief store vectors at the end of the block
 */
#define _fill_store_vectors(_blk) { \
	/* store direction array */ \
	(_blk)->dir = dir; \
	/* store large offset */ \
	(_blk)->offset = offset; \
	/* store diff vectors */ \
	_store(_pv((_blk)->diff.dh), dh); \
	_store(_pv((_blk)->diff.dv), dv); \
	_store(_pv((_blk)->diff.de), de); \
	_store(_pv((_blk)->diff.df), df); \
	/* store delta vectors */ \
	_store(_pv((_blk)->sd.delta), delta); \
	_store(_pv((_blk)->sd.max), max); \
}

/**
 * @fn fill_test_xdrop
 * @brief returns negative if terminate-condition detected
 */
static _force_inline
int64_t fill_test_xdrop(
	struct sea_dp_context_s *this,
	struct sea_block_s *blk)
{
	return(this->tx - blk->sd.max[BW/2]);
}

/**
 * @fn fill_bulk_test_ij_bound
 * @brief returns negative if ij-bound (for the bulk fill) is invaded
 */
static _force_inline
int64_t fill_bulk_test_ij_bound(
	struct sea_dp_context_s *this,
	struct sea_block_s *blk)
{
	return 0;
}

/**
 * @fn fill_cap_test_ij_bound
 * @brief returns negative if ij-bound (for the cap fill) is invaded
 */
static _force_inline
int64_t fill_cap_test_ij_bound(
	struct sea_dp_context_s *this,
	struct sea_block_s *blk)
{
	return 0;
}

/**
 * @fn fill_bulk_test_p_bound
 * @brief returns negative if p-bound (for the bulk fill) is invaded
 */
static _force_inline
int64_t fill_bulk_test_p_bound(
	struct sea_dp_context_s *this,
	uint32_t p)
{
	return(this->rr.s.limp - p);
}

/**
 * @fn fill_create_head
 * @brief create joint_head on the stack to start block extension
 */
static _force_inline
struct sea_block_s *fill_create_head(
	struct sea_dp_context_s *this,
	struct sea_joint_tail_s *prev_tail)
{
	/* init working stack */
	struct sea_joint_head_s *head = (struct sea_joint_head_s *)this->stack_top;
	head->tail = prev_tail;

	/* copy phantom vectors */
	struct sea_block_s *blk = _phantom_block(head + 1);
	_memcpy_blk_aa(
		blk + SEA_BLOCK_PHANTOM_OFFSET,
		_last_block(prev_tail) + SEA_BLOCK_PHANTOM_OFFSET,
		SEA_BLOCK_PHANTOM_SIZE);
	return(blk + 1);
}

/**
 * @fn fill_create_tail
 * @brief create joint_tail at the end of the blocks
 */
static _force_inline
struct sea_joint_tail_s *fill_create_tail(
	struct sea_dp_context_s *this,
	struct sea_joint_tail_s *prev_tail,
	struct sea_block_s *blk,
	uint32_t p)
{
	/* create joint_tail */
	struct sea_joint_tail_s *tail = (struct sea_joint_tail_s *)blk;
	tail->v = prev_tail->v;				/* copy middle deltas */

	/* search max section */
	v32i16_t md = _load_v32i16(prev_tail->v);
	v32i16_t sd = _cvt_v32i8_v32i16(_load(&blk->sd));
	int16_t max = _hmax_v32i16(_add_v32i16(md, sd));
	tail->max = max + blk->offset;

	/* save misc to joint_tail */
	tail->p = p;
	tail->mp = -1;
	tail->psum = p + prev_tail->psum;

	/* */

	/* update dp context */
	this->stack_top = (void *)(tail + 1);
	return(tail);
}

/**
 * @fn fill_bulk_block
 * @brief fill a block
 */
static _force_inline
void fill_bulk_block(
	struct sea_dp_context_s *this,
	struct sea_block_s *blk)
{
	/* load vectors onto registers */
	_fill_load_contexts(blk - 1);

	/**
	 * @macro _fill_block
	 * @brief an element of unrolled fill-in loop
	 */
	#define _fill_block(_direction, _label, _jump_to) { \
		if(_unlikely(!_dir_is_##_direction(dir))) { \
			goto _linear_fill_##_jump_to; \
		} \
		_linear_fill_##_label: _fill_##_direction(); \
		if(--i == 0) { break; } \
	}

	/* update diff vectors */
	int64_t i = BLK;
	struct sea_mask_pair_s *mask_ptr = blk->mask;
	while(1) {					/* 4x unrolled loop */
		_fill_block(down, d1, r1);
		_fill_block(right, r1, d2);
		_fill_block(down, d2, r2);
		_fill_block(right, r2, d1);
	}

	/* update seq offset */
	_fill_update_offset();

	/* store vectors */
	_fill_store_vectors(blk);

	return;
}

/**
 * @struct sea_fill_status_s
 * @brief result container for block fill functions
 */
struct sea_fill_status_s {
	struct sea_block_s *blk;
	uint32_t stat;
	uint32_t p;
};

/**
 * @fn fill_bulk_predetd_blocks
 * @brief fill <blk_cnt> contiguous blocks without ij-bound test
 */
static _force_inline
struct sea_fill_status_s fill_bulk_predetd_blocks(
	struct sea_dp_context_s *this,
	struct sea_block_s *blk,
	uint64_t blk_cnt)
{
	uint32_t stat = CONT;
	uint64_t bc = 0;
	for(bc = 0; bc < blk_cnt; bc++) {
		/* check xdrop termination */
		if(fill_test_xdrop(this, blk - 1) < 0) {
			stat = TERM; break;
		}

		/* fetch sequence */
		rd_bulk_fetch(this);

		/* bulk fill */
		fill_bulk_block(this, blk++);
	}
	return((struct sea_fill_status_s){ blk, stat, (uint32_t)bc * BLK });
}

/**
 * @fn fill_bulk_seq_bounded
 * @brief fill blocks with ij-bound test
 */
static _force_inline
struct sea_fill_status_s fill_bulk_seq_bounded(
	struct sea_dp_context_s *this,
	struct sea_block_s *blk)
{
	uint32_t stat = CONT;

	/* init local coordinate */
	int64_t p = 0;

	/* bulk fill loop */
	while(1) {
		/* check termination */
		if((fill_test_xdrop(this, blk - 1)
		  | fill_bulk_test_ij_bound(this, blk - 1)
		  | fill_bulk_test_p_bound(this, p)) < 0) {
			break;
		}

		/* fetch sequence */
		rd_bulk_fetch(this);
		
		/* bulk fill */
		fill_bulk_block(this, blk++);
		
		/* update p-coordinate */
		p += BLK;
	}
	if(fill_test_xdrop(this, blk) < 0) { stat = TERM; }
	return((struct sea_fill_status_s){ blk, stat, (uint32_t)p });
}

/**
 * @fn fill_cap_seq_bounded
 * @brief fill blocks with cap test
 */
static _force_inline
struct sea_fill_status_s fill_cap_seq_bounded(
	struct sea_dp_context_s *this,
	struct sea_block_s *blk)
{
	uint32_t stat = CONT;

#if 0
	/* check if ij bound is already invaded */
	if(_fill_test_ij_bound(this, blk) < 0) {
		goto _fill_cap_seq_bounded_finish;
	}
#endif
	uint64_t i = 0;
	int64_t p = 0;

	do {
		/* check xdrop termination */
		if(fill_test_xdrop(this, blk - 1) < 0) {
			stat = TERM; goto _fill_cap_seq_bounded_finish;
		}
		/* fetch sequence */
		rd_cap_fetch(this);

		/**
		 * @macro _fill_block_cap
		 * @brief an element of unrolled fill-in loop
		 */
		#define _fill_block_cap() { \
			if(_dir_is_right(dir)) { \
				_fill_right(); \
			} else { \
				_fill_down(); \
			} \
		}

		/* vectors on registers inside this block */ {
			_fill_load_contexts(blk - 1);

			/* update diff vectors */
			struct sea_mask_pair_s *mask_ptr = blk->mask;
			for(i = 0; i < BLK; i++) {
				_fill_block_cap();
				if(fill_cap_test_ij_bound(this, blk) < 0) {
					/* adjust reminders */
					blk->mask[BLK-1] = blk->mask[i];
					_dir_adjust_reminder(dir, i);
					p -= BLK - i - 1;
					break;
				}
			}
			
			/* update seq offset */
			_fill_update_offset();
			
			/* store mask and vectors */
			_fill_store_vectors(blk);
		}
		
		/* update block pointer and p-coordinate */
		blk++; p += BLK;

	} while(i == BLK);

_fill_cap_seq_bounded_finish:;
	return((struct sea_fill_status_s){ blk, stat, (uint32_t)p });
}

/**
 * @fn calc_max_bulk_blocks_mem
 * @brief calculate maximum number of blocks (limited by stack size)
 */
static _force_inline
uint64_t calc_max_bulk_blocks_mem(
	struct sea_dp_context_s *this)
{
	uint64_t const rem = sizeof(struct sea_joint_head_s)
					   + sizeof(struct sea_joint_tail_s)
					   + 3 * sizeof(struct sea_block_s);
	uint64_t mem_size = this->stack_end - this->stack_top;
	return((mem_size - rem) / sizeof(struct sea_block_s) / BLK);
}

/**
 * @fn calc_max_bulk_blocks_seq
 * @brief calculate maximum number of blocks (limited by seq bounds)
 */
static _force_inline
uint64_t calc_max_bulk_blocks_seq(
	struct sea_dp_context_s *this)
{
	uint64_t max_bulk_p = MIN3(
		this->rr.s.body.alen,
		this->rr.s.body.blen,
		this->rr.s.limp);
	return(max_bulk_p / BLK);
}

/**
 * @fn fill_mem_bounded
 * @brief fill <blk_cnt> contiguous blocks without seq bound tests, adding head and tail
 */
static _force_inline
struct sea_chain_status_s fill_mem_bounded(
	struct sea_dp_context_s *this,
	struct sea_joint_tail_s *prev_tail,
	uint64_t blk_cnt)
{
	struct sea_block_s *blk = fill_create_head(this, prev_tail);
	struct sea_fill_status_s stat = fill_bulk_predetd_blocks(this, blk, blk_cnt);
	struct sea_joint_tail_s *tail = fill_create_tail(this, prev_tail, stat.blk, stat.p);
	return((struct sea_chain_status_s){ tail, stat.stat });
}

/**
 * @fn fill_seq_bounded
 * @brief fill blocks with seq bound tests, adding head and tail
 */
static _force_inline
struct sea_chain_status_s fill_seq_bounded(
	struct sea_dp_context_s *this,
	struct sea_joint_tail_s *prev_tail)
{
	struct sea_block_s *blk = fill_create_head(this, prev_tail);
	struct sea_fill_status_s stat = { blk, TERM };

	/* calculate block size */
	uint64_t seq_bulk_blocks = calc_max_bulk_blocks_seq(this);
	while(seq_bulk_blocks > MIN_BULK_BLOCKS) {
		/* bulk fill without ij-bound test */
		stat = fill_bulk_predetd_blocks(this, stat.blk, seq_bulk_blocks);
		if(stat.stat == TERM) {
			goto _fill_seq_bounded_finish;	/* skip cap */
		}
		seq_bulk_blocks = calc_max_bulk_blocks_seq(this);
	}

	/* bulk fill with ij-bound test */
	stat = fill_bulk_seq_bounded(this, blk);
	if(stat.stat == TERM) {
		goto _fill_seq_bounded_finish;	/* skip cap */
	}

	/* cap fill (without p-bound test) */
	stat = fill_cap_seq_bounded(this, stat.blk);

_fill_seq_bounded_finish:;
	struct sea_joint_tail_s *tail = fill_create_tail(this, prev_tail, stat.blk, stat.p);
	return((struct sea_chain_status_s){ tail, stat.stat });
}

/**
 * @fn fill
 * @brief fill dp matrix inside section pairs
 */
int32_t sea_dp_add_stack(
	struct sea_dp_context_s *this);
struct sea_chain_status_s fill(
	struct sea_dp_context_s *this,
	struct sea_joint_tail_s *prev_tail,
	struct sea_section_pair_s *sec)
{
	/* calculate max p-coordinate fit in the remaining stack */
	struct sea_chain_status_s stat = { NULL, TERM };

	/* init section and restore sequence reader buffer */
	rd_load_section(this, sec);
	rd_load_seq(this, prev_tail);

	/* calculate block sizes */
	uint64_t mem_bulk_blocks = calc_max_bulk_blocks_mem(this);
	uint64_t seq_bulk_blocks = calc_max_bulk_blocks_seq(this);

	/* extra large bulk fill (with stack allocation) */
	while(_unlikely(mem_bulk_blocks < seq_bulk_blocks)) {
		if(mem_bulk_blocks > MIN_BULK_BLOCKS) {
			stat = fill_mem_bounded(this, prev_tail, mem_bulk_blocks);
			if(stat.stat == TERM) { goto _fill_finish; }

			/* fill-in area has changed */
			seq_bulk_blocks = calc_max_bulk_blocks_seq(this);
		}

		/* malloc the next stack and set pointers */
		sea_dp_add_stack(this);

		/* stack size has changed */
		mem_bulk_blocks = calc_max_bulk_blocks_mem(this);
	}

	/* bulk fill with seq bound check */
	stat = fill_seq_bounded(this, prev_tail);

_fill_finish:;
	/* restore section */
	rd_save_section(this, sec);
	rd_save_seq(this, (struct sea_joint_tail_s *)stat.ptr);

	return(stat);
}

#if 0
/**
 * @val dynamic, guided
 */
static
struct sea_aln_funcs_s dynamic[2] = {
	{ NULL, NULL, NULL },
	{ wide_dynamic_fill, NULL, NULL }
	// {narrow_dynamic_fill, narrow_dynamic_merge, narrow_dynamic_trace},
	// {wide_dynamic_fill, wide_dynamic_merge, wide_dynamic_trace}
};
static
struct sea_aln_funcs_s guided[2] = {
	{ NULL, NULL, NULL },
	{ NULL, NULL, NULL }
	// {narrow_guided_fill, narrow_guided_merge, narrow_guided_trace},
	// {wide_guided_fill, wide_guided_merge, wide_guided_trace}
};
#endif

/**
 * @fn extract_max
 * @brief extract max value from 8-bit 16-cell vector
 */
static _force_inline
int8_t extract_max(int8_t const vector[][4])
{
	int8_t *v = (int8_t *)vector;
	int8_t max = -128;
	for(int i = 0; i < 16; i++) {
		max = (v[i] > max) ? v[i] : max;
	}
	return(max);
}

/**
 * @fn sea_init_restore_default_params
 */
static _force_inline
void sea_init_restore_default_params(
	struct sea_params_s *params)
{
	#define restore(_name, _default) { \
		params->_name = ((uint64_t)(params->_name) == 0) ? (_default) : (params->_name); \
	}
	// restore(band_width, 		SEA_WIDE);
	// restore(band_type, 			SEA_DYNAMIC);
	restore(seq_a_format, 		SEA_ASCII);
	restore(seq_a_direction, 	SEA_FW_ONLY);
	restore(seq_b_format, 		SEA_ASCII);
	restore(seq_b_direction, 	SEA_FW_ONLY);
	restore(aln_format, 		SEA_ASCII);
	restore(head_margin, 		0);
	restore(tail_margin, 		0);
	restore(xdrop, 				100);
	restore(score_matrix, 		SEA_SCORE_SIMPLE(1, 1, 1, 1));
	return;
}

/**
 * @fn sea_init_create_score_vector
 */
static _force_inline
struct sea_score_vec_s sea_init_create_score_vector(
	struct sea_score_s const *score_matrix)
{
	#define broadcast(x) { \
		(x), (x), (x), (x), \
		(x), (x), (x), (x), \
		(x), (x), (x), (x), \
		(x), (x), (x), (x) \
	}

	int8_t *v = (int8_t *)score_matrix->score_sub;
	struct sea_score_vec_s sc;
	for(int i = 0; i < 16; i++) {
		sc.sbv[i] = v[i];
		sc.geav[i] = -score_matrix->score_ge_a;
		sc.gebv[i] = -score_matrix->score_ge_b;
		sc.giav[i] = -score_matrix->score_gi_a;
		sc.gibv[i] = -score_matrix->score_gi_b;
	}
	return(sc);
}

/**
 * @fn sea_init_create_dir_dynamic
 */
static _force_inline
union sea_dir_u sea_init_create_dir_dynamic(
	struct sea_score_s const *score_matrix)
{
	return((union sea_dir_u) {
		.dynamic = {
			0,				/* zero independent of scoreing schemes */
			0x80000000		/* (0, 0) -> (0, 1) */
		}
	});
}

/**
 * @fn sea_init_create_small_delta
 */
static _force_inline
struct sea_small_delta_s sea_init_create_small_delta(
	struct sea_score_s const *score_matrix)
{
	int8_t max = extract_max(score_matrix->score_sub);
	int8_t diff_a = max + score_matrix->score_ge_a;
	int8_t diff_b = -score_matrix->score_ge_b;

	struct sea_small_delta_s sd;
	for(int i = 0; i < BW/2; i++) {
		sd.delta[i] = diff_a;
		sd.delta[BW/2 + i] = diff_b;
		sd.max[i] = 0;
		sd.max[BW/2 + i] = -diff_b;
	}
	return(sd);
}

/**
 * @fn sea_init_create_middle_delta
 */
static _force_inline
struct sea_middle_delta_s sea_init_create_middle_delta(
	struct sea_score_s const *score_matrix)
{
	int8_t max = extract_max(score_matrix->score_sub);
	int16_t coef_a = -max - 2*score_matrix->score_ge_a;
	int16_t coef_b = -max - 2*score_matrix->score_ge_b;
	int16_t ofs_a = -score_matrix->score_gi_a;
	int16_t ofs_b = -score_matrix->score_gi_b;

	struct sea_middle_delta_s md;
	for(int i = 0; i < BW/2; i++) {
		md.delta[i] = ofs_a + coef_a * (BW/2 - i);
		md.delta[BW/2 + i] = ofs_b + coef_b * i;
	}
	md.delta[BW/2] = 0;
	return(md);
}

/**
 * @fn sea_init_create_diff_vectors
 */
static _force_inline
struct sea_diff_vec_s sea_init_create_diff_vectors(
	struct sea_score_s const *score_matrix)
{
	int8_t max = extract_max(score_matrix->score_sub);
	int8_t drop_dh = 0;
	int8_t raise_dh = max + 2*score_matrix->score_ge_b;
	int8_t drop_dv = 0;
	int8_t raise_dv = max + 2*score_matrix->score_ge_a;
	int8_t drop_de = -score_matrix->score_gi_a + score_matrix->score_ge_a;
	int8_t drop_df = -score_matrix->score_gi_b + score_matrix->score_ge_b;

	struct sea_diff_vec_s diff;
	for(int i = 0; i < BW/2; i++) {
		diff.dh[i] = drop_dh;
		diff.dh[BW/2 + i] = raise_dh;
		diff.dv[i] = raise_dv;
		diff.dv[BW/2 + i] = drop_dv;
		diff.de[i] = drop_de;
		diff.df[i] = drop_df;
	}
	return(diff);
}

/**
 * @fn sea_init
 */
sea_t *sea_init(
	struct sea_params_s const *params)
{
	/* sequence reader table */
	void (*const rd_table[3][7])(
		uint8_t *dst,
		uint8_t const *src,
		uint64_t idx,
		uint64_t src_len,
		uint64_t copy_len) = {
		[SEA_FW_ONLY] = {
			[SEA_ASCII] = _load_ascii_fw,
			// [SEA_4BIT] = _load_4bit_fw,
			// [SEA_2BIT] = _load_2bit_fw,
			// [SEA_4BIT8PACKED] = _load_4bit8packed_fw,
			// [SEA_2BIT8PACKED] = _load_2bit8packed_fw
		},
		[SEA_FW_RV] = {
			[SEA_ASCII] = _load_ascii_fr,
			// [SEA_4BIT] = _load_4bit_fr,
			// [SEA_2BIT] = _load_2bit_fr,
			// [SEA_4BIT8PACKED] = _load_4bit8packed_fr,
			// [SEA_2BIT8PACKED] = _load_2bit8packed_fr
		}
	};
	/* alignment writer table */
	struct sea_writer_s wr_fw_table[4] = {
		[SEA_STR] = {
			.push = _push_ascii_f,
			.type = WR_ASCII,
			.fr = WR_FW
		},
		[SEA_CIGAR] = {
			.push = _push_cigar_f,
			.type = WR_CIGAR,
			.fr = WR_FW
		},
		[SEA_DIR] = {
			.push = _push_dir_f,
			.type = WR_DIR,
			.fr = WR_FW
		}
	};
	struct sea_writer_s wr_rv_table[4] = {
		[SEA_STR] = {
			.push = _push_ascii_r,
			.type = WR_ASCII,
			.fr = WR_RV
		},
		[SEA_CIGAR] = {
			.push = _push_cigar_r,
			.type = WR_CIGAR,
			.fr = WR_RV
		},
		[SEA_DIR] = {
			.push = _push_dir_r,
			.type = WR_DIR,
			.fr = WR_RV
		}
	};

	if(params == NULL) {
		debug("params must not be NULL");
		return(NULL);
	}

	/* copy params to local stack */
	struct sea_params_s params_intl = *params;

	/* restore defaults */
	sea_init_restore_default_params(&params_intl);

	/* malloc sea_context */
	struct sea_context_s *ctx = (struct sea_context_s *)sea_aligned_malloc(
		sizeof(struct sea_context_s),
		SEA_MEM_ALIGN_SIZE);
	if(ctx == NULL) {
		return(NULL);
	}

	/* fill context */
	*ctx = (struct sea_context_s) {
		/* template */
		.k = (struct sea_dp_context_s) {
			.stack_top = NULL,						/* filled on dp init */
			.stack_end = NULL,						/* filled on dp init */
			.pdr = NULL,							/* filled on dp init */
			.tdr = NULL,							/* filled on dp init */

			.ll = (struct sea_writer_work_s) { 0 },	/* work: no need to init */
			.rr = (struct sea_reader_work_s) { 0 },	/* work: no need to init */
			.l = (struct sea_writer_s)(ctx->rv),	/* reverse writer (default on init) */
			.r = (struct sea_reader_s) {
				.loada = rd_table
					[params_intl.seq_a_direction]
					[params_intl.seq_a_format],		/* seq a reader */
				.loadb = rd_table
					[params_intl.seq_b_direction]
					[params_intl.seq_b_format]		/* seq b reader */
			},

			.scv = sea_init_create_score_vector(params_intl.score_matrix),
			.tx = params_intl.xdrop,
			// .max = 0,								/* at (0, 0) */
			// .m_tail = &ctx->tail,					/* at (0, 0) */

			.mem_cnt = 0,
			.mem_size = INIT_STACK_SIZE,
			.mem_array = { 0 },		/* NULL */

			// .fn = ((params_intl.band_type == SEA_GUIDED) ? guided : dynamic)[0],
		},
		.md = sea_init_create_middle_delta(params_intl.score_matrix),
		.blk = (struct sea_phantom_block_s) {
			.mask = {
				{ 0x00000000, 0x00000000 },
				{ 0x0000ffff, 0xffff0000 }
			},
			.dir = sea_init_create_dir_dynamic(params_intl.score_matrix),
			.offset = 0,
			.diff = sea_init_create_diff_vectors(params_intl.score_matrix),
			.sd = sea_init_create_small_delta(params_intl.score_matrix)
		},
		.tail = (struct sea_joint_tail_s) {
			.v = &ctx->md,
			.p = 2,
			.mp = 0,
			.mq = 0,
			.psum = 2,
			.wa = { 0 },
			.wb = { 0 }
		},
		.rv = wr_rv_table[params_intl.aln_format],
		.fw = wr_fw_table[params_intl.aln_format],
		.params = params_intl
	};

	return((sea_t *)ctx);
}

/**
 * @fn sea_clean
 */
void sea_clean(
	struct sea_context_s *ctx)
{
	if(ctx != NULL) { sea_aligned_free(ctx); }
	return;
}

/**
 * @fn sea_dp_init
 */
struct sea_dp_context_s *sea_dp_init(
	struct sea_context_s const *ctx,
	struct sea_seq_pair_s const *p,
	uint8_t const *guide,
	uint64_t glen)
{
	/* malloc stack memory */
	struct sea_dp_context_s *this = (struct sea_dp_context_s *)sea_aligned_malloc(
		ctx->k.mem_size,
		SEA_MEM_ALIGN_SIZE);
	if(this == NULL) {
		debug("failed to malloc memory");
		return(NULL);
	}

	/* init seq pointers */
	_memcpy_blk_aa(
		&this->rr.p,
		p,
		sizeof(struct sea_seq_pair_s));

	/* load template */
	_memcpy_blk_aa(
		(uint8_t *)this + SEA_DP_CONTEXT_LOAD_OFFSET,
		(uint8_t *)&ctx->k + SEA_DP_CONTEXT_LOAD_OFFSET,
		SEA_DP_CONTEXT_LOAD_SIZE);

	/* init stack pointers */
	this->stack_top = (uint8_t *)this + sizeof(struct sea_dp_context_s);
	this->stack_end = (uint8_t *)this + this->mem_size;
	this->pdr = guide;
	this->tdr = guide + glen;
	this->tail = (struct sea_joint_tail_s *)&ctx->tail;

	return(this);
}

/**
 * @fn sea_dp_build_stat
 */
struct sea_chain_status_s sea_dp_build_stat(
	struct sea_dp_context_s const *this)
{
	return((struct sea_chain_status_s){ this->tail, SEA_SUCCESS });
}

/**
 * @fn sea_dp_add_stack
 */
int32_t sea_dp_add_stack(
	struct sea_dp_context_s *this)
{
	uint8_t *ptr = (uint8_t *)sea_aligned_malloc(
		(this->mem_size *= 2),
		SEA_MEM_ALIGN_SIZE);
	if(ptr == NULL) {
		this->mem_size /= 2;
		return(SEA_ERROR_OUT_OF_MEM);
	}
	this->mem_array[this->mem_cnt++] = this->stack_top = ptr;
	this->stack_end = this->stack_top + this->mem_size;
	return(SEA_SUCCESS);
}

/**
 * @fn sea_dp_malloc
 */
void *sea_dp_malloc(
	struct sea_dp_context_s *this,
	uint64_t size)
{
	/* roundup */
	uint64_t const align_size = 16;
	size += (align_size - 1);
	size &= ~(align_size - 1);

	/* malloc */
	if((this->stack_end - this->stack_top) < size) {
		if(this->mem_size < size) { this->mem_size = size; }
		if(sea_dp_add_stack(this) != SEA_SUCCESS) {
			return(NULL);
		}
	}
	this->stack_top += size;
	return((void *)(this->stack_top - size));
}

/**
 * @fn sea_dp_free
 */
void sea_dp_free(
	struct sea_dp_context_s *this,
	void *ptr)
{
	/* nothing to do */
	return;
}

/**
 * @fn sea_dp_clean
 */
void sea_dp_clean(
	struct sea_dp_context_s *this)
{
	if(this == NULL) {
		return;
	}

	for(uint64_t i = 0; i < SEA_MEM_ARRAY_SIZE; i++) {
		sea_aligned_free(this->mem_array[i]);
	}
	sea_aligned_free(this);
	return;
}

/**
 * @fn sea_align_dynamic
 */
struct sea_result *sea_align_dynamic(
	struct sea_context_s const *ctx,
	struct sea_seq_pair_s const *seq,
	struct sea_checkpoint_s const *cp,
	uint64_t cplen)
{
	return(NULL);
}

/**
 * @fn sea_align_guided
 */
struct sea_result *sea_align_guided(
	struct sea_context_s const *ctx,
	struct sea_seq_pair_s const *seq,
	struct sea_checkpoint_s const *cp,
	uint64_t cplen,
	uint8_t const *guide,
	uint64_t glen)
{
	return(NULL);
}

/* unittests */
#ifdef TEST
#include <assert.h>
void unittest(void)
{
	char const *a = "AAAAAAAAAAAAAAAA";
	char const *b = "AAAAAAAAAAAAAAAA";

	/* sea_init */
	sea_t *ctx = sea_init(SEA_PARAMS(
		.seq_a_format = SEA_ASCII,
		.seq_a_direction = SEA_FW_ONLY,
		.seq_b_format = SEA_ASCII,
		.seq_b_direction = SEA_FW_ONLY,
		.aln_format = SEA_ASCII,
		.xdrop = 100,
		.score_matrix = SEA_SCORE_SIMPLE(1, 1, 1, 1)));
	assert(ctx != NULL);

	/* build dp context */
	sea_seq_pair_t seq = sea_build_seq_pair(a, strlen(a), b, strlen(b));

	sea_dp_t *dp = sea_dp_init(ctx, &seq, NULL, 0);
	assert(dp != NULL);

	dump(dp, 1024);

	struct sea_chain_status_s stat = sea_dp_build_stat(dp);
	dump(stat.ptr - sizeof(struct sea_joint_tail_s) - sizeof(struct sea_phantom_block_s) - sizeof(struct sea_joint_head_s), 1024);

	struct sea_section_pair_s sec = sea_build_section_pair(
		sea_build_section(0, 16, 0, 16),
		sea_build_section(16, 16, 16, 16),
		32);

	dump(&sec, 80);

	stat = fill(dp, stat.ptr, &sec);

	dump(stat.ptr, sizeof(struct sea_joint_tail_s));

	sea_dp_clean(dp);

	/* sea_clean */
	sea_clean(ctx);

	return;
}

#endif

#ifdef MAIN
#include <stdio.h>
int main(void)
{
	unittest();
	return 0;
}
#endif

/**
 * end of branch2.c
 */
