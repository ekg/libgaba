/**
 * @file graph_example.c
 */

// #include "gaba.h"			/* single target configuration: gcc graph_example.c gaba.c -DMODEL=AFFINE -DBW=64 */
//#include "gaba_wrap.h"			/* multiple target configuration: gcc graph_example.c libgaba.a */
#include "gaba.h"
#include "gaba_parse.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

int printer(void *fp, int64_t len, char c)
{
	return(fprintf((FILE *)fp, "%ld%c", len, c));
}

char* strto4bit(char* s) {
    uint64_t l = strlen(s);
    char* b = malloc(l*sizeof(char));
    for (uint64_t i = 0; i < l; ++i) {
        char k = '\x00';
        switch (s[i]) {
        case 'A': k = '\x01'; break;
        case 'C': k = '\x02'; break;
        case 'G': k = '\x04'; break;
        case 'T': k = '\x08'; break;
        default: break;
        }
        b[i] = k;
    }
    return b;
}

typedef struct {
    struct gaba_section_s sec;
    struct gaba_section_s query;
    struct gaba_fill_s fill;
    struct gaba_node_s** prev;
    uint16_t count_prev;
    struct gaba_node_s** next;
    uint16_t count_next;
} _gaba_node_s;
typedef _gaba_node_s gaba_node_s;

typedef struct {
    uint32_t size;
    gaba_node_s** nodes;
} gaba_graph_s;

gaba_node_s* gaba_create_node(uint32_t id, const char* seq) {
    gaba_node_s* node = calloc(1, sizeof(gaba_node_s));
    node->sec = gaba_build_section(id, seq, strlen(seq));
    return node;
}

void gaba_node_add_prev(gaba_node_s* n, gaba_node_s* m) {
    ++n->count_prev;
    n->prev = realloc(n->prev, n->count_prev*sizeof(gaba_node_s*));
    n->prev[n->count_prev -1] = m;
}

void gaba_node_add_next(gaba_node_s* n, gaba_node_s* m) {
    ++n->count_prev;
    n->next = (gaba_node_s**)realloc(n->next, n->count_next*sizeof(gaba_node_s*));
    n->next[n->count_next -1] = m;
}

int32_t gaba_graph_add_node(gaba_graph_s* graph, gaba_node_s* node) {
    if (graph->size % 1024 == 0) {
        size_t old_size = graph->size * sizeof(void*);
        size_t increment = 1024 * sizeof(void*);
        if (!(graph->nodes = realloc((void*)graph->nodes, old_size + increment))) {
            fprintf(stderr, "error:[gaba] could not allocate memory for graph\n"); exit(1);
        }
    }
    ++graph->size;
    graph->nodes[graph->size-1] = node;
    return graph->size;
}

void gaba_graph_align(gaba_graph_s* graph, char* q) {
    /* create config */
    gaba_t *ctx = gaba_init(
        GABA_PARAMS(
            .xdrop = 100,
            GABA_SCORE_SIMPLE(1, 4, 6, 1)					/* match award, mismatch penalty, gap open penalty (G_i), and gap extension penalty (G_e) */
            ));
    char const t[64] = { 0 };  /* tail array */
    struct gaba_section_s tail = gaba_build_section(4, t, 64);
    // for each node, fill
    // if we have no predecessors, fill root
    // if we have some, merge them and then fill
    // each node needs to maintain its own section and the query section (awkward?)
	/* create thread-local object */
	void const *lim = (void const *)0x800000000000;		/* end-of-userland pointer */
	gaba_dp_t *dp = gaba_dp_init(ctx, lim, lim);		/* dp[0] holds a 64-cell-wide context */


    
    struct gaba_fill_s const *f = gaba_dp_fill_root(dp,	/* dp -> &dp[_dp_ctx_index(band_width)] makes the band width selectable */
                                                    ap, 0,  /* a-side (reference side) sequence and start position */
                                                    bp, 0,	/* b-side (query) */
                                                    UINT32_MAX);  /* max extension length */


	struct gaba_fill_s const *m = f;					/* track max */
    
	while((f->status & GABA_TERM) == 0) {
		if(f->status & GABA_UPDATE_A) { ap = &tail; }	/* substitute the pointer by the tail section's if it reached the end */
		if(f->status & GABA_UPDATE_B) { bp = &tail; }

		f = gaba_dp_fill(dp, f, ap, bp, UINT32_MAX);	/* extend the banded matrix */
		m = f->max > m->max ? f : m;					/* swap if maximum score was updated */
	}

	struct gaba_alignment_s *r = gaba_dp_trace(dp,
		m,												/* section with the max */
		NULL											/* custom allocator: see struct gaba_alloc_s in gaba.h */
	);

	printf("score(%ld), path length(%lu)\n", r->score, r->plen);
	gaba_print_cigar_forward(
		printer, (void *)stdout,						/* printer */
		r->path,										/* bit-encoded path array */
		0,												/* offset is always zero */
		r->plen											/* path length */
	);
	printf("\n");

    
}

int main(int argc, char *argv[]) {

	/* init section pointers */
	struct gaba_section_s const *ap = &asec, *bp = &bsec;
	struct gaba_fill_s const *f = gaba_dp_fill_root(dp,	/* dp -> &dp[_dp_ctx_index(band_width)] makes the band width selectable */
		ap, 0,											/* a-side (reference side) sequence and start position */
		bp, 0,											/* b-side (query) */
		UINT32_MAX										/* max extension length */
	);

	/* until X-drop condition is detected */
	struct gaba_fill_s const *m = f;					/* track max */
	while((f->stat & GABA_TERM) == 0) {
		if(f->stat & GABA_UPDATE_A) { ap = &tail; }		/* substitute the pointer by the tail section's if it reached the end */
		if(f->stat & GABA_UPDATE_B) { bp = &tail; }

		f = gaba_dp_fill(dp, f, ap, bp, UINT32_MAX);	/* extend the banded matrix */
		m = f->max > m->max ? f : m;					/* swap if maximum score was updated */
	}

	struct gaba_alignment_s *r = gaba_dp_trace(dp,
		m,												/* section with the max */
		NULL											/* custom allocator: see struct gaba_alloc_s in gaba.h */
	);

	printf("score(%ld), path length(%lu)\n", r->score, r->plen);
	gaba_dp_print_cigar_forward(
		printer, (void *)stdout,						/* printer */
		r->path,										/* bit-encoded path array */
		0,												/* offset is always zero */
		r->plen											/* path length */
	);
	printf("\n");

	/* clean up */
	gaba_dp_res_free(r);
	gaba_dp_clean(dp);
	gaba_clean(ctx);
    free(a); free(b);
	return 0;
}

/**
 * end of example.c
 */
