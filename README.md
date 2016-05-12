# libgaba

GABA provides affine-gap penalty adaptive semi-global banded alignment on nucleotide ~~string graphs~~ (current implementation only supports trees). It uses fixed-width adaptive banded semi-global (a variant of Smith-Waterman and Gotoh) algorithm combined with difference DP (an acceleration technique similar to the Myers' bit-parallel editdist algorithm) and X-drop termination.

## Usage

```
/* init context */
gaba_t *ctx = gaba_init(GABA_PARAMS(
	.xdrop = 100,
	.score_matrix = GABA_SCORE_SIMPLE(2, 3, 5, 1)));

/* create section info (gaba_build_section is a macro) */
struct gaba_section_s asec = gaba_build_section(1, 0, strlen(a));
struct gaba_section_s bsec = gaba_build_section(2, 0, strlen(b));

/*
 * lim points the end of memory region of forward
 * sequences. If the reverse sequences are not provided,
 * lim should be 0x800000000000 (the tail address
 * of the user space)
 */
void const *lim = (void const *)0x800000000000;

/* init dp context */
gaba_dp_t *dp = gaba_dp_init(ctx, lim, lim);

/* fill root with asec and bsec from (0, 0) */
struct gaba_section_s *ap = &asec, *bp = &bsec;
struct gaba_fill_s *f = gaba_dp_fill_root(dp, ap, 0, bp, 0);

/*
 * f->status & GABA_STATUS_UPDATE_A indicates section a
 * reached the end.
 */
if(f->status & GABA_STATUS_UPDATE_A) {
	ap = /* pointer to the next section of a */
}
if(f->status & GABA_STATUS_UPDATE_B) {
	bp = /* pointer to the next section of b */
}

/* sections are filled from head to tail */
f = gaba_dp_fill(dp, fill, ap, bp);

/*
 * Each fill object (f) contains the max score of the
 * section (f->max). All the fill objects are allocated
 * from the dp context so you do not have to (must not)
 * free them even though you want to discard them.
 */

/*
 * Traces (alignment paths) are calculated from
 * arbitrary sections to the root. Traces are stored
 * in r->path in the array of 1-bit direction
 * (0: RIGHT / 1: DOWN). The generated path and the
 * reversed path of the third argument (reverse section)
 * are concatenated at the root.
 */
struct gaba_result_s *r = gaba_dp_trace(dp,
	f,			/* forward section */
	NULL,		/* reverse section */
	NULL);

/*
 * 1-bit direction string can be converted to CIGAR
 * string (defined in the SAM format) with dump_cigar
 * or print_cigar function.
 */
gaba_dp_print_cigar(
	(gaba_dp_fprintf_t)fprintf,
	(void *)stdout,
	r->path->array,
	r->path->offset,
	r->path->len);

/* destroy the dp context */
gaba_dp_clean(dp);

/* destroy the global context */
gaba_clean(ctx);
```


## Functions

### Init / cleanup global context

#### gaba\_init

Global context holds global configurations of the alignment. It can be shared between multiple threads (thread-safe) and must be created once at first.

```
gaba_t *gaba_init(gaba_params_t const *params);
```

#### gaba\_clean

```
void gaba_clean(gaba_t *ctx);
```

### Init / cleanup DP context

#### gaba\_dp\_init

DP context is a local context holding thread-local configurations and working buffers. It can't be shared between threads.

```
struct gaba_dp_context_s *gaba_dp_init(
	gaba_t const *ctx,
	gaba_seq_pair_t const *p);
```

#### gaba\_dp\_flush

Flush the working buffer in the DP context. Any previous result (generated by `gaba_dp_fill` and `gaba_dp_trace`) will be invalid after call.

```
void gaba_dp_flush(
	gaba_dp_t *this,
	gaba_seq_pair_t const *p);
```

#### gaba\_dp\_clean

```
void gaba_dp_clean(
	gaba_dp_t *this);
```

### Alignment functions

#### gaba\_dp\_fill\_root

Create a root section of the extension alignment.

```
gaba_fill_t *gaba_dp_fill_root(
	gaba_dp_t *this,
	gaba_section_t const *a,
	uint32_t apos,
	gaba_section_t const *b,
	uint32_t bpos);
```

#### gaba\_dp\_fill

Extends a section.

```
gaba_fill_t *gaba_dp_fill(
	gaba_dp_t *this,
	gaba_fill_t const *prev_sec,
	gaba_section_t const *a,
	gaba_section_t const *b);
```

#### gaba\_dp\_merge

Merge multiple alignment sections. (not implemented yet)

```
gaba_fill_t *gaba_dp_merge(
	gaba_dp_t *this,
	gaba_fill_t const *sec_list,
	uint64_t tail_list_len);
```

#### gaba\_dp\_trace

Traceback from the max. score in the given sections and concatenate paths at the root.

```
gaba_result_t *gaba_dp_trace(
	gaba_dp_t *this,
	gaba_fill_t const *fw_tail,
	gaba_fill_t const *rv_tail,
	gaba_clip_params_t const *clip);
```

### Utils

#### gaba\_dp\_print\_cigar

Convert path string and dump to `FILE *`.

```
typedef int (*gaba_dp_fprintf_t)(void *, char const *, ...);
int64_t gaba_dp_print_cigar(
	gaba_dp_fprintf_t fprintf,
	void *fp,
	uint32_t const *path,
	uint32_t offset,
	uint32_t len);
```

#### gaba\_dp\_dump\_cigar

Convert path string and dump to a `buf`. `buf` must have enough room to store the cigar string.

```
int64_t gaba_dp_dump_cigar(
	char *buf,
	uint64_t buf_size,
	uint32_t const *path,
	uint32_t offset,
	uint32_t len);
```

## License

Apache v2.

Copyright (c) 2016 Hajime Suzuki