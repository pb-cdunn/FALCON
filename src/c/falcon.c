/*
 *
 *       Filename:  fastcon.c
 *
 *    Description:
 *
 *        Version:  0.1
 *        Created:  07/20/2013 17:00:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jason Chin,
 *        Company:
 *

 #################################################################################$$
 # Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
 #
 # All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted (subject to the limitations in the
 # disclaimer below) provided that the following conditions are met:
 #
 #  * Redistributions of source code must retain the above copyright
 #  notice, this list of conditions and the following disclaimer.
 #
 #  * Redistributions in binary form must reproduce the above
 #  copyright notice, this list of conditions and the following
 #  disclaimer in the documentation and/or other materials provided
 #  with the distribution.
 #
 #  * Neither the name of Pacific Biosciences nor the names of its
 #  contributors may be used to endorse or promote products derived
 #  from this software without specific prior written permission.
 #
 # NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
 # GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
 # BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 # WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 # OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 # DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 # LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 # USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 # OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 # OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 # SUCH DAMAGE.
 #################################################################################$$
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "common.h"

// #define DEBUG_DETAILED_VERBOSE
// #define DEBUG_PRINT_CONS_STATUS

typedef struct {
    seq_coor_t t_pos;
    uint8_t delta;
    char q_base;
    seq_coor_t p_t_pos;   // the tag position of the previous base
    uint8_t p_delta; // the tag delta of the previous base
    char p_q_base;        // the previous base
    unsigned q_id;
} align_tag_t;

typedef struct {
    seq_coor_t len;
    align_tag_t * align_tags;
} align_tags_t;


typedef struct {
    uint16_t size;
    uint16_t n_link;
    seq_coor_t * p_t_pos;   // the tag position of the previous base
    uint8_t * p_delta; // the tag delta of the previous base
    char * p_q_base;        // the previous base
    uint16_t * link_count;
    uint16_t count;
    seq_coor_t best_p_t_pos;
    uint8_t best_p_delta;
    uint8_t best_p_q_base; // encoded base
    double score;
} align_tag_col_t;

typedef struct {
    align_tag_col_t * base;
} msa_base_group_t;

typedef struct {
    uint8_t size;
    uint8_t max_delta;
    msa_base_group_t * delta;
} msa_delta_group_t;

typedef msa_delta_group_t * msa_pos_t;

align_tags_t * get_align_tags( char * aln_q_seq,
                               char * aln_t_seq,
                               seq_coor_t aln_seq_len,
                               seq_coor_t q_start,
                               seq_coor_t t_start,
                               unsigned q_id,
                               seq_coor_t t_offset) {
    char p_q_base;
    align_tags_t * tags;
    seq_coor_t i, j, jj, k, p_j, p_jj;

    tags = calloc( 1, sizeof(align_tags_t) );
    tags->len = aln_seq_len;
    tags->align_tags = calloc( aln_seq_len + 1, sizeof(align_tag_t) );
    i = q_start - 1;
    j = t_start - 1;
    jj = 0;
    p_j = -1;
    p_jj = 0;
    p_q_base = '.';

    for (k = 0; k < aln_seq_len; k++) {
        if (aln_q_seq[k] != '-') {
            i ++;
            jj ++;
        }
        if (aln_t_seq[k] != '-') {
            j ++;
            jj = 0;
        }
        //printf("t %d %d %d %c %c\n", q_id, j, jj, aln_t_seq[k], aln_q_seq[k]);


        if ( j + t_offset >= 0 && jj < UINT8_MAX && p_jj < UINT8_MAX) {
            (tags->align_tags[k]).t_pos = j + t_offset;
            (tags->align_tags[k]).delta = jj;
            (tags->align_tags[k]).p_t_pos = p_j + t_offset;
            (tags->align_tags[k]).p_delta = p_jj;
            (tags->align_tags[k]).p_q_base = p_q_base;
            (tags->align_tags[k]).q_base = aln_q_seq[k];
            (tags->align_tags[k]).q_id = q_id;

            p_j = j;
            p_jj = jj;
            p_q_base = aln_q_seq[k];
        } else {
            // If this happens, that's fine. We simply break on large insertions,
            // and the tag array will be ended below the loop at this k.
            break; // when there is a big alignment gap > UINT8_MAX, stop extending the tagging string
        }
    }
    // sentinal at the end
    //k = aln_seq_len;
    tags->len = k;
    (tags->align_tags[k]).t_pos = UINT_MAX;
    (tags->align_tags[k]).delta = UINT8_MAX;
    (tags->align_tags[k]).q_base = '.';
    (tags->align_tags[k]).q_id = UINT_MAX;
    return tags;
}

void free_align_tags( align_tags_t * tags) {
    free( tags->align_tags );
    free( tags );
}


void allocate_aln_col( align_tag_col_t * col) {
    col->p_t_pos = ( seq_coor_t * ) calloc(col->size, sizeof( seq_coor_t ));
    col->p_delta = ( uint8_t * ) calloc(col->size, sizeof( uint8_t ));
    col->p_q_base = ( char * )calloc(col->size, sizeof( char ));
    col->link_count = ( uint16_t * ) calloc(col->size, sizeof( uint16_t ));
}

void realloc_aln_col( align_tag_col_t * col ) {
    col->p_t_pos = (seq_coor_t *) realloc( col->p_t_pos, (col->size) * sizeof( seq_coor_t ));
    col->p_delta = ( uint8_t *)  realloc( col->p_delta, (col->size) * sizeof( uint8_t ));
    col->p_q_base = (char *) realloc( col->p_q_base, (col->size) * sizeof( char ));
    col->link_count = ( uint16_t *) realloc( col->link_count, (col->size) * sizeof( uint16_t ));
}

void free_aln_col( align_tag_col_t * col) {
    free(col->p_t_pos);
    free(col->p_delta);
    free(col->p_q_base);
    free(col->link_count);
}


void allocate_delta_group( msa_delta_group_t * g) {
    int i,j;
    g->max_delta = 0;
    g->delta = (msa_base_group_t *) calloc( g->size, sizeof(msa_base_group_t));
    for (i = 0; i< g->size; i++) {
        g->delta[i].base = ( align_tag_col_t * ) calloc( 5, sizeof(align_tag_col_t ) );
        for (j = 0; j < 5; j++ ) {
             g->delta[i].base[j].size = 8;
             allocate_aln_col(&(g->delta[i].base[j]));
        }
    }
}

void realloc_delta_group( msa_delta_group_t * g, uint16_t new_size ) {
    int i, j, bs, es;
    bs = g->size;
    es = new_size;
    g->delta = (msa_base_group_t *) realloc(g->delta, new_size * sizeof(msa_base_group_t));
    for (i=bs; i < es; i++) {
        g->delta[i].base = ( align_tag_col_t *) calloc( 5, sizeof(align_tag_col_t ) );
        for (j = 0; j < 5; j++ ) {
             g->delta[i].base[j].size = 8;
             allocate_aln_col(&(g->delta[i].base[j]));
        }
    }
    g->size = new_size;
}

void free_delta_group( msa_delta_group_t * g) {
    //manything to do here
    int i, j;
    for (i = 0; i < g->size; i++) {
        for (j = 0; j < 5; j++) {
            free_aln_col( &(g->delta[i].base[j]) );
        }
        free(g->delta[i].base);
    }
    free(g->delta);
}

void update_col( align_tag_col_t * col, seq_coor_t p_t_pos, uint8_t p_delta, char p_q_base) {
    int updated = 0;
    int kk;
    col->count += 1;
    for (kk = 0; kk < col->n_link; kk++) {
        if ( p_t_pos == col->p_t_pos[kk] &&
             p_delta == col->p_delta[kk] &&
             p_q_base == col->p_q_base[kk] ) {
            col->link_count[kk] ++;
            updated = 1;
            break;
        }
    }
    if (updated == 0) {
        if (col->n_link + 1 > col->size) {
            if (col->size < (UINT16_MAX >> 1)-1) {
                col->size *= 2;
            } else {
                col->size += 256;
            }
            if (col->size >= UINT16_MAX-1) {
                fprintf(stderr, "[update_col] Assert will fail! (col->size < (UINT16_MAX-1))? Values: %u >= %u\n", col->size, (UINT16_MAX-1));
            }
            // fflush(stderr);
            assert( col->size < UINT16_MAX-1 );
            realloc_aln_col(col);
        }
        kk = col->n_link;

        col->p_t_pos[kk] = p_t_pos;
        col->p_delta[kk] = p_delta;
        col->p_q_base[kk] = p_q_base;
        col->link_count[kk] = 1;
        col->n_link++;
    }
}


msa_pos_t * get_msa_working_sapce(unsigned int max_t_len) {
    msa_pos_t * msa_array;
    unsigned int i;
    msa_array = calloc(max_t_len, sizeof(msa_pos_t *));
    for (i = 0; i < max_t_len; i++) {
        msa_array[i] = calloc(1, sizeof(msa_delta_group_t));
        msa_array[i]->size = 8;
        allocate_delta_group(msa_array[i]);
    }
    return msa_array;
}

void clean_msa_working_space( msa_pos_t * msa_array, unsigned int max_t_len) {
    unsigned int i,j,k;
    align_tag_col_t * col;
    for (i = 0; i < max_t_len; i++) {
        for (j =0; j < msa_array[i]->max_delta + 1; j++) {
            for (k = 0; k < 5; k++ ) {
                col = msa_array[i]->delta[j].base + k;
                /*
                for (c =0; c < col->size; c++) {
                    col->p_t_pos[c] = 0;
                    col->p_delta[c] = 0;
                    col->p_q_base[c] = 0;
                    col->link_count[c] =0;
                }
                */
                col->n_link = 0;
                col->count = 0;
                col->best_p_t_pos = 0;
                col->best_p_delta = 0;
                col->best_p_q_base = 0;
                col->score = 0;
            }
        }
        msa_array[i]->max_delta = 0;
    }
}

#define STATIC_ALLOCATE
//#undef STATIC_ALLOCATE

consensus_data * get_cns_from_align_tags( align_tags_t ** tag_seqs,
                                          unsigned n_tag_seqs,
                                          unsigned t_len,
                                          unsigned min_cov ) {

    seq_coor_t i, j;
    seq_coor_t t_pos = 0;
    unsigned int * coverage;
    unsigned int * local_nbase;

    consensus_data * consensus;
    //char * consensus;
    align_tag_t * c_tag;

    coverage = calloc( t_len, sizeof(unsigned int) );
    local_nbase = calloc( t_len, sizeof(unsigned int) );

#ifndef STATIC_ALLOCATE

    msa_pos_t * msa_array = NULL; // For more efficiency, this should be injected.
    msa_array = calloc(t_len, sizeof(msa_pos_t *));

    for (i = 0; i < t_len; i++) {
        msa_array[i] = calloc(1, sizeof(msa_delta_group_t));
        msa_array[i]->size = 8;
        allocate_delta_group(msa_array[i]);
    }

#else

    unsigned const max_t_len = 128000;
    if (t_len > max_t_len) {
        fprintf(stderr, "t_len==%d > %d\n", t_len, max_t_len);
        //abort();
        return 0;
    }
    static msa_pos_t * msa_array = NULL;
    if ( msa_array == NULL) {
        msa_array = get_msa_working_sapce(max_t_len + 1);
    }

    if (t_len > max_t_len) {
        fprintf(stderr, "[get_cns_from_align_tags]  Assert will fail!\n");
    }
    assert(t_len <= max_t_len);

#endif

    // loop through every alignment
    //printf("XX %d\n", n_tag_seqs);
    for (i = 0; i < n_tag_seqs; i++) {

        // for each alignment position, insert the alignment tag to msa_array
        for (j = 0; j < tag_seqs[i]->len; j++) {
            c_tag = tag_seqs[i]->align_tags + j;
            unsigned int delta;
            delta = c_tag->delta;
            if (delta == 0) {
                t_pos = c_tag->t_pos;
                coverage[ t_pos ] ++;
            }
            // Assume t_pos was set on earlier iteration.
            // (Otherwise, use its initial value, which might be an error. ~cd)
            if (delta > msa_array[t_pos]->max_delta) {
                msa_array[t_pos]->max_delta = delta;
                if (msa_array[t_pos]->max_delta + 4 > msa_array[t_pos]->size ) {
                    realloc_delta_group(msa_array[t_pos], msa_array[t_pos]->max_delta + 8);
                }
            }

            unsigned int base = -1;
            switch (c_tag->q_base) {
                case 'A': base = 0; break;
                case 'C': base = 1; break;
                case 'G': base = 2; break;
                case 'T': base = 3; break;
                case '-': base = 4; break;
                default:
                    base = -1;
                    fprintf(stderr, "WARNING: Bad input detected! c_tag->q_base = '%c' (int value = %d).\n", c_tag->q_base, (int) c_tag->q_base);
                    fprintf(stderr, "[get_cns_from_align_tags]:   before update_col: j = %d / %d, i = %d\n", j, tag_seqs[i]->len, i);
                    fprintf(stderr, "t_pos = %d, delta = %d, base = %d, c_tag->q_base = %c, c_tag->p_t_pos = %d, c_tag->p_delta = %d, c_tag->p_q_base = %d\n",
                                t_pos, delta, base, c_tag->q_base, c_tag->p_t_pos, c_tag->p_delta, c_tag->p_q_base);
                    break;
            }

            // Note: On bad input, base may be -1.
            update_col( &(msa_array[t_pos]->delta[delta].base[base]), c_tag->p_t_pos, c_tag->p_delta, c_tag->p_q_base);
            local_nbase[ t_pos ] ++;
        }
    }

#ifdef DEBUG_DETAILED_VERBOSE
    fprintf(stderr, "[get_cns_from_align_tags] 3: Ping!\n");
#endif

    // propogate score throught the alignment links, setup backtracking information
    align_tag_col_t * g_best_aln_col = 0;
    unsigned int g_best_k = 0;
    seq_coor_t g_best_t_pos = 0;
    {
        int k;
	int best_k;
        double score;
        double best_score;
        double g_best_score;

        align_tag_col_t * aln_col;

        g_best_score = -1;

        for (i = 0; i < t_len; i++) {  //loop through every template base
            #ifdef DEBUG_DETAILED_VERBOSE
                fprintf(stderr, "[get_cns_from_align_tags] 3.1a: starting i = %d / %d\n", i, t_len);
            #endif

            //printf("max delta: %d %d\n", i, msa_array[i]->max_delta);
            for (j = 0; j <= msa_array[i]->max_delta; j++) { // loop through every delta position
                #ifdef DEBUG_DETAILED_VERBOSE
                    fprintf(stderr, "     j = %d / %d; k = ", j, msa_array[i]->max_delta);
                #endif
                for (k = 0; k < 5; k++) {  // loop through diff bases of the same delta posiiton
                    aln_col = msa_array[i]->delta[j].base + k;
                    #ifdef DEBUG_DETAILED_VERBOSE
                        fprintf(stderr, "%d ", k);
                    #endif
                    if (aln_col->count >= 0) {
                        best_score = -1;

                        for (int link = 0; link < aln_col->n_link; link++) { // loop through differnt link to previous column
                            int pi;
                            int pj;
                            int pk;
                            pi = aln_col->p_t_pos[link];
                            pj = aln_col->p_delta[link];
                            switch (aln_col->p_q_base[link]) {
                                case 'A': pk = 0; break;
                                case 'C': pk = 1; break;
                                case 'G': pk = 2; break;
                                case 'T': pk = 3; break;
                                case '-': pk = 4; break;
                                default: pk = 4;
                            }

                            if (aln_col->p_t_pos[link] == -1) {
                                score =  (double) aln_col->link_count[link] - (double) coverage[i] * 0.5;
                            } else {
                                score = msa_array[pi]->delta[pj].base[pk].score +
                                        (double) aln_col->link_count[link] - (double) coverage[i] * 0.5;
                            }
                            if (score > best_score) {
                                best_score = score;
                                aln_col->best_p_t_pos = pi;
                                aln_col->best_p_delta = pj;
                                aln_col->best_p_q_base = pk;
				best_k = k;
                            }
                        }
                        aln_col->score = best_score;
                        if (best_score > g_best_score) {
                            g_best_score = best_score;
                            g_best_aln_col = aln_col;
                            g_best_t_pos = i;
			    g_best_k = best_k;
                        }
                    }
                }
                #ifdef DEBUG_DETAILED_VERBOSE
                    fprintf(stderr, "\n");
                #endif
            }

            #ifdef DEBUG_DETAILED_VERBOSE
                fprintf(stderr, "[get_cns_from_align_tags] 3.1b: ending i = %d\n", i);
            #endif
        }
        if (g_best_score == -1) {
            fprintf(stderr, "In get_cns_from_align_tags(), g_best_score==-1\n");
            return 0;
        }

        if (g_best_score == -1) {
            // This can't actually happen, because there is a check right above.
            fprintf(stderr, "[get_cns_from_align_tags]  Assert will fail!\n");
        }
        assert(g_best_score > 0);
    }

    #ifdef DEBUG_DETAILED_VERBOSE
        fprintf(stderr, "[get_cns_from_align_tags] 4: Ping!\n");
    #endif

    // reconstruct the sequences
    unsigned int index;
    char bb = '$';
    int k;
    char * cns_str;
    int * eqv;
    double score0;

    consensus = calloc( 1, sizeof(consensus_data) );
    consensus->sequence = calloc( t_len * 2 + 1, sizeof(char) );
    consensus->eqv = calloc( t_len * 2 + 1, sizeof(unsigned int) );
    cns_str = consensus->sequence;
    eqv =  consensus->eqv;

    index = 0;
    k = g_best_k;
    i = g_best_t_pos;

    while (1) {
        if (coverage[i] > min_cov) {
            switch (k) {
                case 0: bb = 'A'; break;
                case 1: bb = 'C'; break;
                case 2: bb = 'G'; break;
                case 3: bb = 'T'; break;
                case 4: bb = '-'; break;
            }
        } else {
            switch (k) {
                case 0: bb = 'a'; break;
                case 1: bb = 'c'; break;
                case 2: bb = 'g'; break;
                case 3: bb = 't'; break;
                case 4: bb = '-'; break;
            }
        }
        // Note: On bad input, bb will keep previous value, possibly '$'.

        score0 = g_best_aln_col->score;
        i = g_best_aln_col->best_p_t_pos;
        if (i == -1 || index >= t_len * 2) break;
        j = g_best_aln_col->best_p_delta;
        k = g_best_aln_col->best_p_q_base;
        g_best_aln_col = msa_array[i]->delta[j].base + k;

        if (bb != '-') {
            cns_str[index] = bb;
            eqv[index] = (int) score0 - (int) g_best_aln_col->score;
            //printf("C %d %d %c %lf %d %d\n", i, index, bb, g_best_aln_col->score, coverage[i], eqv[index] );
            index ++;
        }
    }

    #ifdef DEBUG_DETAILED_VERBOSE
        fprintf(stderr, "[get_cns_from_align_tags] 5: Ping!\n");
    #endif

    // reverse the sequence
    for (i = 0; i < index/2; i++) {
        cns_str[i] = cns_str[i] ^ cns_str[index-i-1];
        cns_str[index-i-1] = cns_str[i] ^ cns_str[index-i-1];
        cns_str[i] = cns_str[i] ^ cns_str[index-i-1];
        eqv[i] = eqv[i] ^ eqv[index-i-1];
        eqv[index-i-1] = eqv[i] ^ eqv[index-i-1];
        eqv[i] = eqv[i] ^ eqv[index-i-1];
    }

    cns_str[index] = 0;
    //printf("%s\n", cns_str);
#ifndef STATIC_ALLOCATE
    for (i = 0; i < t_len; i++) {
        free_delta_group(msa_array[i]);
        free(msa_array[i]);
    }

    free(msa_array);
#else
    clean_msa_working_space(msa_array, t_len+1);
#endif

    free(coverage);
    free(local_nbase);

    #ifdef DEBUG_DETAILED_VERBOSE
        fprintf(stderr, "[get_cns_from_align_tags] 6: Ping!\n");
    #endif

    return consensus;
}

//const unsigned int K = 8;

consensus_data * generate_consensus( char ** input_seq,
                           unsigned int n_seq,
                           unsigned min_cov,
                           unsigned K,
                           double min_idt) {
    unsigned int j;
    unsigned int seq_count;
    unsigned int aligned_seq_count;
    kmer_lookup * lk_ptr;
    seq_array sa_ptr;
    seq_addr_array sda_ptr;
    kmer_match * kmer_match_ptr;
    aln_range * arange;
    alignment * aln;
    align_tags_t ** tags_list;
    //char * consensus;
    consensus_data * consensus;
    double max_diff;
    const unsigned int lk_ptr_size = (1 << (K * 2));
    max_diff = 1.0 - min_idt;

    fprintf(stderr, "[consensus] In generate_consensus.\n");

    seq_count = n_seq;
    //printf("XX n_seq %d\n", n_seq);
    //for (j=0; j < seq_count; j++) {
    //    printf("seq_len: %u %u\n", j, strlen(input_seq[j]));
    //};
    fflush(stdout);

    tags_list = calloc( seq_count, sizeof(align_tags_t *) );
    lk_ptr = allocate_kmer_lookup( lk_ptr_size );
    sa_ptr = allocate_seq( (seq_coor_t) strlen( input_seq[0]) );
    sda_ptr = allocate_seq_addr( (seq_coor_t) strlen( input_seq[0]) );
    add_sequence( 0, K, input_seq[0], strlen(input_seq[0]), sda_ptr, sa_ptr, lk_ptr);
    mask_k_mer(lk_ptr_size, lk_ptr, 10000);

    aligned_seq_count = 0;
    for (j=1; j < seq_count; j++) {

        //printf("seq_len: %ld %u\n", j, strlen(input_seq[j]));

        kmer_match_ptr = find_kmer_pos_for_seq(input_seq[j], strlen(input_seq[j]), K, sda_ptr, lk_ptr);
#define INDEL_ALLOWENCE_0 6

        arange = find_best_aln_range(kmer_match_ptr, K, K * INDEL_ALLOWENCE_0, 5);  // narrow band to avoid aligning through big indels

        //printf("1:%ld %ld %ld %ld\n", arange_->s1, arange_->e1, arange_->s2, arange_->e2);

        //arange = find_best_aln_range2(kmer_match_ptr, K, K * INDEL_ALLOWENCE_0, 5);  // narrow band to avoid aligning through big indels

#ifdef DEBUG_PRINT_CONS_STATUS
        fprintf(stderr, "(internal) 1:%ld %ld %ld %ld\n", arange->s1, arange->e1, arange->s2, arange->e2);
#endif

#define INDEL_ALLOWENCE_1 0.10
        if (arange->e1 - arange->s1 < 100 || arange->e2 - arange->s2 < 100 ||
            abs( (arange->e1 - arange->s1 ) - (arange->e2 - arange->s2) ) >
                   (int) (0.5 * INDEL_ALLOWENCE_1 * (arange->e1 - arange->s1 + arange->e2 - arange->s2))) {
            free_kmer_match( kmer_match_ptr);
            free_aln_range(arange);
            continue;
        }
        //printf("%ld %s\n", strlen(input_seq[j]), input_seq[j]);
        //printf("%ld %s\n\n", strlen(input_seq[0]), input_seq[0]);


#define INDEL_ALLOWENCE_2 150

        aln = align(input_seq[j]+arange->s1, arange->e1 - arange->s1 ,
                    input_seq[0]+arange->s2, arange->e2 - arange->s2 ,
                    INDEL_ALLOWENCE_2, 1);

#ifdef DEBUG_PRINT_CONS_STATUS
        fprintf(stderr, "(internal) 2: %lf\n\n", (((double) aln->dist / (double) aln->aln_str_size)));
#endif

        if (aln->aln_str_size > 500 && ((double) aln->dist / (double) aln->aln_str_size) < max_diff) {
            tags_list[aligned_seq_count] = get_align_tags( aln->q_aln_str,
                                                           aln->t_aln_str,
                                                           aln->aln_str_size,
                                                           arange->s1,
                                                           arange->s2,
                                                           j,
                                                           0);
            aligned_seq_count ++;
        }
        /***
        for (k = 0; k < tags_list[j]->len; k++) {
            printf("%ld %d %c\n", tags_list[j]->align_tags[k].t_pos,
                                   tags_list[j]->align_tags[k].delta,
                                   tags_list[j]->align_tags[k].q_base);
        }
        ***/
        free_aln_range(arange);
        free_alignment(aln);
        free_kmer_match( kmer_match_ptr);
    }

    if (aligned_seq_count > 0) {
        consensus = get_cns_from_align_tags( tags_list, aligned_seq_count, strlen(input_seq[0]), min_cov );
        if (!consensus) return 0;
    } else {
        // allocate an empty consensus sequence
        consensus = calloc( 1, sizeof(consensus_data) );
        consensus->sequence = calloc( 1, sizeof(char) );
        consensus->eqv = calloc( 1, sizeof(unsigned int) );
    }
    //free(consensus);
    free_seq_addr_array(sda_ptr);
    free_seq_array(sa_ptr);
    free_kmer_lookup(lk_ptr);
    for (j=0; j < aligned_seq_count; j++) {
        free_align_tags(tags_list[j]);
    }
    free(tags_list);
    return consensus;
}

consensus_data * generate_consensus_from_mapping( char ** input_seq,
                           aln_range **input_aranges,
                           unsigned int n_seq,
                           unsigned min_cov,
                           unsigned K,
                           double min_idt) {
    unsigned int j = 0;
    unsigned int aln_clip_offset = 0;
    unsigned int seq_count = 0;
    unsigned int aligned_seq_count = 0;
    kmer_match * kmer_match_ptr = NULL;
    aln_range * arange = NULL;
    alignment * aln = NULL;
    align_tags_t ** tags_list = NULL;
    //char * consensus;
    consensus_data * consensus = NULL;
    double max_diff = 0.0;
    max_diff = 1.0 - min_idt;

    fprintf(stderr, "[consensus] In generate_consensus_from_mapping.\n");

    seq_count = n_seq;

    fflush(stdout);

    tags_list = calloc( seq_count, sizeof(align_tags_t *) );

    aligned_seq_count = 0;
    for (j=1; j < seq_count; j++) {
        arange = input_aranges[j];

#define INDEL_ALLOWENCE_1 0.10

#ifdef DEBUG_PRINT_CONS_STATUS
        fprintf(stderr, "(external) 1: j = %d / %d, %ld %ld %ld %ld\n", j, seq_count, arange->s1, arange->e1, arange->s2, arange->e2);
#endif

        if (arange->e1 - arange->s1 < 100 || arange->e2 - arange->s2 < 100 ||
            abs( (arange->e1 - arange->s1 ) - (arange->e2 - arange->s2) ) >
                   (int) (0.5 * INDEL_ALLOWENCE_1 * (arange->e1 - arange->s1 + arange->e2 - arange->s2))) {
            continue;
        }

#define INDEL_ALLOWENCE_2 150
        aln = align(input_seq[j]+arange->s1, arange->e1 - arange->s1 ,
                    input_seq[0]+arange->s2, arange->e2 - arange->s2 ,
                    INDEL_ALLOWENCE_2, 1);

#ifdef DEBUG_PRINT_CONS_STATUS
        fprintf(stderr, "(external) 2: %lf\n", (((double) aln->dist / (double) aln->aln_str_size)));
#endif

#ifdef DEBUG_PRINT_ALN
        fprintf(stderr, "T: %s\n", aln->t_aln_str);
        fprintf(stderr, "Q: %s\n", aln->q_aln_str);
#endif

        // Find the first non-indel aligned base and offset the coordinates.
        seq_coor_t q_start = arange->s1;
        seq_coor_t t_start = arange->s2;
        for (aln_clip_offset = 0; aln_clip_offset < aln->aln_str_size; ++aln_clip_offset) {
            if (aln->t_aln_str[aln_clip_offset] != '-' && aln->q_aln_str[aln_clip_offset] != '-') {
                break;
            }
            if (aln->q_aln_str[aln_clip_offset] != '-') {
                ++q_start;
            }
            if (aln->t_aln_str[aln_clip_offset] != '-') {
                ++t_start;
            }
        }

#ifdef DEBUG_PRINT_CONS_STATUS
        if (aln_clip_offset > 0) {
            fprintf(stderr, "(external) 2.1: aln_clip_offset = %u, aln->aln_str_size = %d\n", aln_clip_offset, aln->aln_str_size);
        }
#endif

        // if (aln->t_aln_str[] == '-' || aln->q_aln_str[0] == '-' || aln->t_aln_str[aln->aln_str_size-1] == '-' || aln->q_aln_str[aln->aln_str_size-1] == '-') {
        //     fprintf(stderr, "(external) 2.1: Skipping adding of alignment because it has leading/trailing insertions.\n");
        //     continue;
        // }

        if ((aln->aln_str_size - aln_clip_offset) > 500 && ((double) aln->dist / (double) (aln->aln_str_size - aln_clip_offset)) < max_diff) {
            tags_list[aligned_seq_count] = get_align_tags( aln->q_aln_str + aln_clip_offset,
                                                           aln->t_aln_str + aln_clip_offset,
                                                           aln->aln_str_size - aln_clip_offset,
                                                           q_start,
                                                           t_start,
                                                           j,
                                                           0);
            aligned_seq_count ++;
#ifdef DEBUG_PRINT_CONS_STATUS
            fprintf(stderr, "(external) 3: Passed filters and added tags.\n");
#endif
        }

#ifdef DEBUG_PRINT_CONS_STATUS
        fprintf(stderr, "(external) 4: aligned_seq_count = %d\n\n", aligned_seq_count);
#endif

        free_alignment(aln);
    }

#ifdef DEBUG_PRINT_CONS_STATUS
    fprintf(stderr, "(external) 5: Finally, aligned_seq_count = %d\n\n", aligned_seq_count);
#endif

    if (aligned_seq_count > 0) {
        consensus = get_cns_from_align_tags( tags_list, aligned_seq_count, strlen(input_seq[0]), min_cov );
        if (!consensus) {
            return 0;
        }
    } else {
        // allocate an empty consensus sequence
        consensus = calloc( 1, sizeof(consensus_data) );
        consensus->sequence = calloc( 1, sizeof(char) );
        consensus->eqv = calloc( 1, sizeof(unsigned int) );
    }

    for (j=0; j < aligned_seq_count; j++) {
        free_align_tags(tags_list[j]);
    }
    free(tags_list);
    return consensus;
}

void free_consensus_data( consensus_data * consensus ){
    free(consensus->sequence);
    free(consensus->eqv);
    free(consensus);
}
