
/**
 * @file guided.h
 *
 * @brief direction determiner and address calculation macros for the guided band algorithm
 */
#ifndef _GUIDED_H_INCLUDED
#define _GUIDED_H_INCLUDED

#include "../sea.h"
#include "../util/util.h"
#include "../arch/dir.h"
#include <stdint.h>

/**
 * @macro BLK
 * @brief block split length
 */
#define BLK 		( 32 )

/**
 * address calculation macros
 */
#define guided_dr_size()		( 0 )
#define guided_blk_num(p, q) ( \
	((p) & ~(BLK-1)) / BLK \
)
#define guided_blk_addr(p, q) ( \
	((p) & (BLK-1)) * bpl() + ((q)+BW/2) * sizeof(cell_t) \
)
#define guided_addr(p, q) ( \
	  guided_blk_num(p, q) * bpb() \
	+ guided_blk_addr(p, q) \
	+ head_size() \
)

/**
 * direction determiner variants
 */
struct _dir {
	uint8_t *pdr;
	uint8_t d2;
};
typedef struct _dir dir_t;

/**
 * @macro dir, dir2, dir_ue, dir_le, dir2_ue, dir2_le
 * @brief direction (or 2-bit direction) of the upper (lower) edge
 */
#define guided_dir_ue(r)	( 0x04 & (r).d2 )
#define guided_dir2_ue(r) 	( 0x05 & (r).d2 )
#define guided_dir_le(r)	( (0x08 & (r).d2)>>1 )
#define guided_dir2_le(r) 	( (0x0a & (r).d2)>>1 )
#define guided_dir(r)		dir_ue(r)
#define guided_dir2(r) 		dir2_ue(r)
#define guided_dir_raw(r) 	( (r).d2 )

#if 0
#define guided_dir_topq(r)		( (int64_t)-1 + (0x01 & ((r).d2>>2)) )
#define guided_dir_leftq(r)		( (int64_t)0x01 & ((r).d2>>2) )
#define guided_dir_topleftq(r)	( (int64_t)(0x01 & (r).d2) - (0x01 & ((~(r).d2)>>2)) )
#endif

/**
 * direction determiners for fill-in macros
 */
/**
 * @macro guided_dir_init
 */
#define guided_dir_init(r, k, dp, p) { \
	(r).pdr = (k)->pdr; \
	(r).d2 = ((r).pdr[p]<<2) | (r).pdr[p-1]; \
	/*(r).d2 = (r).pdr[(p)-1]<<2;*/	/** SEA_LEFT or SEA_TOP */ \
}

/**
 * @macro guided_dir_start_block
 */
#define guided_dir_start_block(r, k, dp, p) { \
	/** nothing to do */ \
}

/**
 * @macro guided_dir_det_next
 * @brief set 2-bit direction flag from direction array.
 */
#define guided_dir_det_next(r, k, dp, p) { \
	uint8_t d = (r).pdr[++(p)]; \
	(r).d2 = (d<<2) | ((r).d2>>2); \
}

/**
 * @macro guided_dir_empty
 */
#define guided_dir_empty(r, k, dp, p) { \
	/** nothing to do */ \
}

/**
 * @macro guided_dir_end_block
 */
#define guided_dir_end_block(r, k, dp, p) { \
	/** nothing to do */ \
}

/**
 * @macro guided_dir_test_bound
 */
#define guided_dir_test_bound(r, k, dp, p) ( \
	((k)->tdr - (r).pdr + 2) - p - BLK \
)
#define guided_dir_test_bound_cap(r, k, dp, p) ( \
	((k)->tdr - (r).pdr + 2) - p \
)

/**
 * direction loaders for search and traceback functions
 */
/**
 * @macro guided_dir_set_pdr
 */
#define guided_dir_set_pdr(r, k, dp, p, sp) { \
	(r).pdr = (k)->pdr; \
	(r).d2 = ((r).pdr[p]<<2) | (r).pdr[p-1]; \
}

/**
 * @macro guided_dir_load_forward
 */
#define guided_dir_load_forward(r, k, dp, p, sp) { \
	(r).d2 = ((r).pdr[++(p)]<<2) | ((r).d2>>2); \
}
#define guided_dir_go_forward(r, k, dp, p, sp) { \
	p++; \
}

/**
 * @macro guided_dir_load_backward
 * @brief set 2-bit direction flag in reverse access.
 */
#define guided_dir_load_backward(r, k, dp, p, sp) { \
	p--; \
	(r).d2 = 0x0f & (((r).d2<<2) | (r).pdr[p - 1]); \
}
#define guided_dir_go_backward(r, k, dp, p, sp) { \
	p--; \
}

/**
 * @macro guided_dir_sum_i_blk
 * @brief calculate sum of diff_i from p to the end of the block
 */
#define guided_dir_sum_i_blk(r, k, dp, p, sp) ( \
	dir_vec_sum_i((r).pdr + (sp) + (((p) - (sp)) & ~(BLK-1)), \
		((p) - (sp)) & (BLK-1)) \
)

/**
 * fast direction loaders
 */
#define guided_dir_set_pdr_fast(r, k, dp, p, sp) { \
	guided_dir_set_pdr(r, k, dp, p, sp); \
}
#define guided_dir_load_backward_fast(r, k, dp, p, sp) { \
	guided_dir_load_backward(r, k, dp, p, sp); \
}

#endif /* #ifndef _GUIDED_H_INCLUDED */
/**
 * end of guided.h
 */