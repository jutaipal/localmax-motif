#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <inttypes.h>
#include <string.h>
#include <time.h>

/* VERSION 0.2, 13 Jan 2025   part of AUTOSEED code modified to enable bitstream input (sequences of arbitrary length)   */
/* bitstream and multithreaded code written with assistance from Claude 3.5 Sonnet and ChatGPT4 and o1                   */
/* Debugged manually by J Taipale as bit shifting is not their forte                                                     */
/* Compile with clang -march=native -ffast-math -O3 -o localmax-motif localmax-motif.c                             */

/* GLOBAL VARIABLES */
uint64_t mask_ULL[42][42];   /* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
uint64_t lowmask_ULL[42];    /* LOW MASK FOR EACH KMER */
uint64_t highmask_ULL[42];   /* HIGH MASK FOR EACH KMER */
short int Nlength = 32;              /* LENGTH OF MASK */

#define NUM_FILES 2
#define MAX_GAP_LENGTH 10
#define KMER_VARIATION 1
#define BUFFER_SIZE 1024
#define MAX_KMER_LEN 32  // Adjust if needed for your use case
#define MAX_SEQ_LEN 500  // Adjust if needed for your use case
#define MAX_WIDTH_OF_PWM 40
#define FLANK_WIDTH 5

struct match {short int *position; double *score;};
short int match_init (struct match *i, short int width)
{
short int maximum_width = Nlength+10;
short int counter;
(*i).position = malloc(sizeof(short int) * maximum_width + 5);
(*i).score = malloc(sizeof(double) * maximum_width + 5);
for (counter = 0; counter < maximum_width; counter++)
{
(*i).position[counter] = 0;
(*i).score[counter] = 0;
}
return(0);
}

/* COUNT PWM */
struct count_pwm {char *name; short int width; long int max_counts; double **incidence;};
short int count_pwm_clear (struct count_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = MAX_WIDTH_OF_PWM;
short int counter;
short int counter2;
strcpy ((*i).name, name);
(*i).width = width;
(*i).max_counts = initial_value;
for (counter = 0; counter < 5; counter++)
{
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).incidence[counter][counter2] = initial_value;
}
return(0);
}
short int count_pwm_init (struct count_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = MAX_WIDTH_OF_PWM;
short int counter;
short int counter2;
(*i).name = malloc(1000);
strcpy ((*i).name, name);
(*i).width = width;
(*i).max_counts = initial_value;
(*i).incidence = malloc(sizeof(double *) * 5 + 5);
for (counter = 0; counter < 5; counter++)
{
(*i).incidence[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).incidence[counter][counter2] = initial_value;
}
return(0);
}
short int count_pwm_free (struct count_pwm *i)
{
    short int counter;
    free((*i).name);
    for (counter = 0; counter < 5; counter++) free((*i).incidence[counter]);
    free((*i).incidence);
    return(0);
}

/* NORMALIZED PWM */
struct normalized_pwm {char *name; char *seed; short int width; long int max_counts; double *information_content; short int *original_position; double *position_score; long int *total_counts_for_column; double **fraction; short int negative_values_allowed;};
short int normalized_pwm_init (struct normalized_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = MAX_WIDTH_OF_PWM;
short int counter;
short int counter2;
(*i).negative_values_allowed = 0;
(*i).name = malloc(100);
strcpy ((*i).name, name);
(*i).seed = malloc(1000);
strcpy ((*i).seed, "UNKNOWN");
(*i).width = width;
(*i).max_counts = initial_value;
(*i).fraction = malloc(sizeof(double *) * 5 + 5);
(*i).information_content = malloc(sizeof(double) * maximum_width + 5);
(*i).position_score = malloc(sizeof(double) * maximum_width + 5);
(*i).original_position = malloc(sizeof(short int) * maximum_width + 5);
(*i).total_counts_for_column = malloc(sizeof(long int) * maximum_width + 5);

for (counter = 0; counter < 5; counter++)
{
(*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).fraction[counter][counter2] = initial_value;
}
for (counter2 = 0; counter2 < maximum_width; counter2++)
{
(*i).information_content[counter2] = 0;
(*i).position_score[counter2] = 0;
(*i).original_position[counter2] = counter2;
(*i).total_counts_for_column[counter2] = 0;
}
return(0);
}
short int normalized_pwm_free (struct normalized_pwm *i)
{
short int counter;
free((*i).name);
free((*i).information_content);
free((*i).position_score);
free((*i).total_counts_for_column);
for (counter = 0; counter < 5; counter++) free((*i).fraction[counter]);
free((*i).fraction);
return(0);
}

/* SUBROUTINE THAT GENERATES PWM FROM IUPAC */
short int Iupac_to_pwm(struct normalized_pwm *n, char *searchstring)
{
short int counter;
short int pwm_position;
short int current_match_position;
short int nucleotide_value;
(*n).width =  strlen(searchstring);

char **canbe;
canbe = malloc(sizeof(char *) * 4 + 5);
for (counter = 0; counter < 4; counter++) canbe[counter] = malloc(200);
strcpy (canbe[0], "100010101001111");
strcpy (canbe[1], "010001100110111");
strcpy (canbe[2], "001010010111011");
strcpy (canbe[3], "000101011011101");

char *conversion_string;
conversion_string = malloc(sizeof(char)*3 + 5);
conversion_string[0] = '\0';
conversion_string[1] = '\0';

char *nucleotide_iupac = "ACGTRYMKWSBDHVN";
short int iupac_length = 15;
    
/* BUILDS PWM */
for(pwm_position = 0; pwm_position < (*n).width; pwm_position++)
{
//printf("\n");
for(current_match_position = 0; (current_match_position < iupac_length) & (searchstring[pwm_position] != nucleotide_iupac[current_match_position]); current_match_position++)
;
    //if(current_match_position == iupac_length) {strcat(seed_story, "\n** SEED ERROR: defective IUPAC"); printf("%s", seed_story); fflush(stdout); exit(1);}
for (nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++)
{
conversion_string[0] = canbe[nucleotide_value][current_match_position];
(*n).fraction[nucleotide_value][pwm_position] = atof(conversion_string) - 1;
}
}

for (counter = 0; counter < 4; counter++) free(canbe[counter]);
free(conversion_string);
return(0);
}

// Print count motif
short int PrintMotif(struct count_pwm *multinomial_motif, short int multinomial, int seed_length) {
    //printf("\nDEBUG: count_pwm structure: width=%d\n", multinomial_motif->width);
    // printf("\nMultinomial %i motif counts, attempting to print all positions:\n", multinomial);
    // printf("Debug: PWM address in print: %p\n", (void*)multinomial_motif);
    // Print positions header
    printf("\n");
    for (int i = 1; i < multinomial_motif->width; i++) {
        printf("\t%d", i);
    }
    printf("\n");
    
    // Print each base as a row
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (int base = 0; base < 4; base++) {
        //printf("%c  ", bases[base]);
        for (int i = 1; i < multinomial_motif->width; i++) {
            printf("\t%.0f", multinomial_motif->incidence[base][i]);
        }
        printf("\n");
    }
    return 0;
}

// Print count motif
short int PrintNormalizedMotif(struct normalized_pwm *motif) {
    // Print positions header
    printf("\n   ");
    for (int i = 0; i < motif->width; i++) {
        printf("\t%d", i);
    }
    printf("\n");
    
    // Print each base as a row
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (int base = 0; base < 4; base++) {
        printf("%c  ", bases[base]);
        for (int i = 0; i < motif->width; i++) {
            printf("\t%.0f", motif->fraction[base][i]);
        }
        printf("\n");
    }
    return 0;
}

short int Findpwmmatch_from_bitstream(struct normalized_pwm *p, double cut_off,
                                    uint64_t *bitstream, size_t sequence_length,
                                    struct match *match) {
   
    //PrintNormalizedMotif(p);
    //printf("\nCutoff: %f", cut_off);
    //printf("\nSequence length: %li", sequence_length);
    //printf("\nSeed motif width: %hi", (*p).width);
    
    uint64_t main_reg = 0;
    uint64_t source_reg = 0;
    const int shift_amount = 2;
    short int number_of_matches = 0;
    short int nucleotide;
    size_t current_position = 0;
    double score;
    signed short int pwm_position;
    size_t last_source_block_moved = 0;
    size_t current_block = 0;
    
    size_t bitstream_length = (sequence_length * 2 + 63) / 64;
    uint64_t *current_bitstream = bitstream;

    // Load initial blocks
    if (bitstream_length > 0) {
        main_reg = *current_bitstream++;
        bitstream_length--;
    }
    if (bitstream_length > 0) {
        source_reg = *current_bitstream++;
        bitstream_length--;
    }

    size_t total_blocks = (sequence_length + 31) / 32;
    size_t total_number_of_blocks_left_to_move_to_source = total_blocks - 2;
    size_t remaining_bits_in_main = 64;
    
    // Process all possible positions
    while (current_position + (*p).width <= sequence_length) {
        // Calculate score for current position
        for (score = 0, pwm_position = (*p).width-1; pwm_position >= 0; pwm_position--) {
            nucleotide = (main_reg & mask_ULL[1][pwm_position]) >> (2 * pwm_position);
            score += (*p).fraction[nucleotide][(*p).width-1-pwm_position];
        }
        //printf("\nScore %f at position: %lu", score, sequence_length - current_position - (*p).width);
        // Record match if score exceeds cutoff
        if (score >= cut_off) {
            number_of_matches++;
            (*match).score[number_of_matches] = score;
            (*match).position[number_of_matches] = sequence_length - current_position - (*p).width;
            //printf("\nMatch!! Score %f at position: %hi",score, (*match).position[number_of_matches]);
        }
        
        // Shift bits
        uint64_t new_bits = (source_reg & 0x3) << (64 - shift_amount);
        main_reg = (main_reg >> shift_amount) | new_bits;
        source_reg >>= shift_amount;
        remaining_bits_in_main -= 2;
        current_position++;

        // Handle last block transition
        if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
            total_number_of_blocks_left_to_move_to_source == 0) {
            remaining_bits_in_main = (sequence_length % 32) * 2;
            last_source_block_moved = 1;
        }
        
        // Load next block if needed
        if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
            total_number_of_blocks_left_to_move_to_source > 0) {
            source_reg = *current_bitstream++;
            total_number_of_blocks_left_to_move_to_source--;
            remaining_bits_in_main = 64;
            current_block++;
        }
    }
    
    (*match).position[0] = number_of_matches;
    return number_of_matches;
}

short int Multinomial_add_to_pwm_from_bitstream(struct count_pwm *p, struct normalized_pwm *qp,
                                              short int match_position, double score, double cut_off,
                                              uint64_t *bitstream, uint64_t *temp_bitstream, size_t sequence_length,
                                              struct normalized_pwm *background_pwm) {
    
    //printf("\nDebugging add_to_pwm: seq_length=%zu, match_position=%d, pwm width=%d\n", sequence_length, match_position, (*p).width);
    
    short int pwm_position;
    short int seed_position;
    short int background_position;
    short int query_position = 0;
    short int in_seed;
    uint64_t main_reg = 0;
    uint64_t source_reg = 0;
    const int shift_amount = 2;
    uint64_t nucleotide;
    size_t last_source_block_moved = 0;
    size_t current_block = 0;

    __atomic_add_fetch(&((*p).max_counts), 1, __ATOMIC_SEQ_CST);
    //(*p).max_counts++;
    
    size_t bitstream_length = (sequence_length * 2 + 63) / 64;
    for(int i = 0; i < bitstream_length; i++) temp_bitstream[i] = bitstream[i]; // copies bitstream
    uint64_t *current_bitstream = temp_bitstream;
    
    // Load initial blocks
    if (bitstream_length > 0) {
        main_reg = *current_bitstream++;
        bitstream_length--;
    }
    if (bitstream_length > 0) {
        source_reg = *current_bitstream++;
        bitstream_length--;
    }

    size_t total_blocks = (sequence_length + 31) / 32;  // How many blocks we should have
    size_t total_number_of_blocks_left_to_move_to_source = total_blocks - 2;
    size_t remaining_bits_in_main = 64;
    
    signed short int shifted_match_position = sequence_length - match_position - (*qp).width - 1;
    short int offset = 0;
    if (sequence_length - match_position <= FLANK_WIDTH + (*qp).width) offset = (*qp).width + FLANK_WIDTH - (sequence_length - match_position);
    
    pwm_position = 2 * FLANK_WIDTH + (*qp).width - offset;
    //printf("\nOffset: %i", offset);
    
    for (; shifted_match_position + (*qp).width + FLANK_WIDTH >= 0 && sequence_length > 1; shifted_match_position--, sequence_length--)
    {

        // shift collect bases from -10 until +10 bp flanking the seed match
        if (shifted_match_position < FLANK_WIDTH)
        {
            nucleotide = main_reg & 0x3;
            seed_position = pwm_position - FLANK_WIDTH;
            if (seed_position >= 0 && seed_position < (*qp).width) in_seed = 1;
            else {in_seed = 0; seed_position = 0;} // to prevent segfault
            
            // EXCLUDES CONSENSUS IF THERE IS MISMATCH
            if (in_seed == 1 && score < cut_off + 1 && (*qp).fraction[nucleotide][seed_position] == 0)
            {
                // printf("\nExcluding nucleotide %c from PWM position %i, shifted match %i, seed position %i", "ACGT"[nucleotide], pwm_position, shifted_match_position, seed_position);
                // (*p).incidence[nucleotide][pwm_position]++;
                pwm_position--;
            }
            else
            {
                // printf("\nAdding nucleotide %c to PWM position %i, shifted match %i", "ACGT"[nucleotide], pwm_position, shifted_match_position);
                __atomic_add_fetch(&((*p).incidence[nucleotide][pwm_position]), 1, __ATOMIC_SEQ_CST);
                //(*p).incidence[nucleotide][pwm_position]++;
                pwm_position--;
            }
            }
        
            // Shift bits
            uint64_t new_bits = (source_reg & 0x3) << (64 - shift_amount);
            main_reg = (main_reg >> shift_amount) | new_bits;
            source_reg >>= shift_amount;
            remaining_bits_in_main -= 2;

            // Handle last block transition
            if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
                total_number_of_blocks_left_to_move_to_source == 0) {
                remaining_bits_in_main = (sequence_length % 32) * 2;
                last_source_block_moved = 1;
            }
            
            // Load next block if needed
            if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
                total_number_of_blocks_left_to_move_to_source > 0) {
                source_reg = *current_bitstream++;
                total_number_of_blocks_left_to_move_to_source--;
                remaining_bits_in_main = 64;
                current_block++;
            }

        
        }
    
    return(0);
}

void print_match_sequences(struct match *matches, uint64_t *bitstream, size_t sequence_length) {
    const size_t PRINT_LENGTH = 40;
    char sequence[PRINT_LENGTH + 1];  // +1 for null terminator
    sequence[PRINT_LENGTH] = '\0';
    
    uint64_t main_reg = 0;
    uint64_t source_reg = 0;
    const int shift_amount = 2;
    size_t last_source_block_moved = 0;
    
    // For each match
    for (int i = 1; i <= matches->position[0]; i++) {
        size_t match_pos = matches->position[i];
        size_t start_pos = 0;
        
        // Fill with 'x' first
        memset(sequence, 'x', PRINT_LENGTH);
        
        // Calculate actual sequence start and length
        if (match_pos > PRINT_LENGTH/2) {
            start_pos = match_pos - PRINT_LENGTH/2;
        }
        
        // Initialize bitstream reading
        size_t bitstream_length = (sequence_length * 2 + 63) / 64;
        uint64_t *current_bitstream = bitstream;
        
        // Load initial blocks
        if (bitstream_length > 0) {
            main_reg = *current_bitstream++;
            bitstream_length--;
        }
        if (bitstream_length > 0) {
            source_reg = *current_bitstream++;
            bitstream_length--;
        }

        size_t total_blocks = (sequence_length + 31) / 32;
        size_t total_number_of_blocks_left_to_move_to_source = total_blocks - 2;
        size_t remaining_bits_in_main = 64;
        
        // Skip to start position
        for (size_t pos = 0; pos < start_pos; pos++) {
            uint64_t new_bits = (source_reg & 0x3) << (64 - shift_amount);
            main_reg = (main_reg >> shift_amount) | new_bits;
            source_reg >>= shift_amount;
            remaining_bits_in_main -= 2;

            if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
                total_number_of_blocks_left_to_move_to_source == 0) {
                remaining_bits_in_main = (sequence_length % 32) * 2;
                last_source_block_moved = 1;
            }
            
            if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
                total_number_of_blocks_left_to_move_to_source > 0) {
                source_reg = *current_bitstream++;
                total_number_of_blocks_left_to_move_to_source--;
                remaining_bits_in_main = 64;
            }
        }
        
        // Extract sequence
        for (size_t pos = start_pos; pos < start_pos + PRINT_LENGTH && pos < sequence_length; pos++) {
            uint64_t nucleotide = main_reg & 0x3;
            sequence[pos - start_pos] = "ACGT"[nucleotide];
            
            // Shift to next position if not at end
            if (pos + 1 < start_pos + PRINT_LENGTH && pos + 1 < sequence_length) {
                uint64_t new_bits = (source_reg & 0x3) << (64 - shift_amount);
                main_reg = (main_reg >> shift_amount) | new_bits;
                source_reg >>= shift_amount;
                remaining_bits_in_main -= 2;

                if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
                    total_number_of_blocks_left_to_move_to_source == 0) {
                    remaining_bits_in_main = (sequence_length % 32) * 2;
                    last_source_block_moved = 1;
                }
                
                if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
                    total_number_of_blocks_left_to_move_to_source > 0) {
                    source_reg = *current_bitstream++;
                    total_number_of_blocks_left_to_move_to_source--;
                    remaining_bits_in_main = 64;
                }
            }
        }
        
        printf("%s\n", sequence);
    }
}


/* SUBROUTINE THAT PRINTS KMER (USES NO MASK) */
void Kmerprint (__uint128_t print_sequence_value, short int kmer_length, short int gap_position, short int gap_length)
{
    for (kmer_length--; kmer_length >= 0; kmer_length--, gap_position--)
    {   if (gap_position == 0) for(int i = 0; i < gap_length; i++) printf("n");
        printf("%c","ACGT"[(print_sequence_value >> (kmer_length * 2)) & 3]);
    }
}

/* SUBROUTINE THAT DETERMINES IF GAP IS AT EITHER OF THE CENTER POSITIONS, IF count_also is set to != 1 returns true */
short int Centergap (short int count_also_spaced_kmers, short int kmer_length, short int gap_position)
{
if (count_also_spaced_kmers != 1) return (1);
if (gap_position == kmer_length / 2) return (1);
if (kmer_length % 2 == 1 && gap_position == kmer_length / 2 - 1) return (1);
else return (0);
}

/* SUBROUTINE THAT DETERMINES IF A GAPPED KMER IS A LOCAL MAXIMUM WITHIN HUDDINGE DISTANCE OF 1 */
/* SEE NITTA ET AL. eLIFE 2015 Methods and Supplementary Figure 1 for algorithm description     */
/* https://doi.org/10.7554/eLife.04837.004                                                      */
/* Bases encoded as bits, A = 00, C = 01, G = 10, T = 11                                        */

short int Localmax(long int *****results, short int file_number, short int current_kmer_length, short int shortest_kmer, short int longest_kmer_counted, short int current_gap_position, short int current_gap_length, long int current_kmer, double kmer_length_difference_cutoff)
{
    short int count_also_spaced_kmers = 1;
    short int too_long_kmer = longest_kmer_counted + 1;

    long int kmer1_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
    long int kmer2_incidence;
    long int compared_kmer = current_kmer;
    signed short int position;
    short int counter;
    short int first_half;
    long int lowbit = 1;  // "01" used to toggle low bit of base encoding (between A-C or G-T), shifted to correct base position in loops
    long int highbit = 2; // "10" used to toggle high bit of base encoding (between A-G or C-T), shifted to correct base position in loops
    short int shift;
    short int true_gap_position = current_gap_position;
    short int true_gap_length = current_gap_length;
    short int start;
    short int end;
    short int left = 0;
    short int right = 1;
    signed short int position2;
    
    /* Substitution; HAMMING OF 1, RETURNS 0 IF ANY KMER WITHIN HAMMING OF 1 HAS HIGHER COUNT */
    for(position=0; position < current_kmer_length; position++, lowbit <<= 2, highbit <<= 2)
    {
        compared_kmer = lowbit ^ current_kmer;  // SHIFTS between A-C or G-T of original kmer at position
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);  // IF SUBSTITUTED KMER HAS HIGHER COUNTS, RETURNS 0 TO INDICATE THAT THE QUERY KMER IS NOT LOCAL MAXIMUM
        compared_kmer = highbit ^ current_kmer; // SHIFTS between A-G or C-T of original kmer at position
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer; // SHIFTS between A-C or G-T of previously shifted kmer to cover last of possible 3 substitutions (total effect is A-T, C-G)
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    /* Shift; FULL SHIFT BY ONE */
    shift = (current_gap_length != 0);
    current_gap_position += shift;
    
    /* ONLY LOOK AT FULL SHIFT FOR UNGAPPED KMERS, OR FOR GAPPED KMERS IF GAP CAN SHIFT ALSO (ALL GAP POSITIONS HAVE BEEN COUNTED) */
    if (current_gap_length == 0 || (count_also_spaced_kmers == 2 && current_gap_position < current_kmer_length) || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2 + 1))
    {
    /* COMPARED FULL SHIFT RIGHT (same shift in gap if any), RETURNS 0 IF ANY KMER WITHIN HAMMING OF 1 HAS HIGHER COUNT  */
    compared_kmer = (current_kmer >> 2) & lowmask_ULL[current_kmer_length-1];
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    lowbit >>= 2; highbit >>= 2;
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    current_gap_position -= shift;
    current_gap_position -= shift;
    if (current_gap_length == 0 || (count_also_spaced_kmers == 2 && current_gap_position > 0) || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2))
    {
    /* COMPARED FULL SHIFT LEFT (same shift in gap if any) */
    compared_kmer = (current_kmer << 2) & lowmask_ULL[current_kmer_length-1];
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    lowbit = 1; highbit = 2;
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    current_gap_position += shift;
    lowbit = 1; highbit = 2;
    if (count_also_spaced_kmers != 0)
    {
    if(current_gap_position == 0) current_gap_position = current_kmer_length / 2;
    
    /* Longer Gap; COMPARE TO KMER WITH LONGER GAP */
    current_gap_length++;
    if (current_gap_length < Nlength - current_kmer_length && current_gap_length < MAX_GAP_LENGTH && true_gap_length != 0)
    {
    compared_kmer = (current_kmer & highmask_ULL[current_kmer_length-current_gap_position]) | ((current_kmer << 2) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    
    /* LOOP TO ANALYZE SHIFT OF EITHER HALF, STARTS WITH SECOND HALF */
    for(first_half = 0; first_half < 2; first_half++)
    {
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    
    /* SWITCHES TO FIRST HALF */
    compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length-current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
    }
    }
        
    /* Shorter Gap; COMPARE TO KMER WITH SHORTER GAP */
    current_gap_length--;
    current_gap_length--;

    lowbit = 1; highbit = 2;
    lowbit <<= ((current_kmer_length - true_gap_position - 1)*2); highbit <<= ((current_kmer_length - true_gap_position - 1)*2);
    if (current_gap_length == 0) current_gap_position = 0;
    if (current_gap_length >= 0)
    {
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - true_gap_position]) | ((current_kmer >> 2) & (lowmask_ULL[current_kmer_length-true_gap_position-1]));
        /* LOOP TO ANALYZE SHIFT OF EITHER HALF, STARTS WITH SECOND HALF */
        for(first_half = 0; first_half < 2; first_half++)
        {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            
            /* SWITCHES TO FIRST HALF */
            compared_kmer = lowmask_ULL[current_kmer_length-1] & ((current_kmer << 2) & highmask_ULL[current_kmer_length - true_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-true_gap_position-1]));
            lowbit <<= 2; highbit <<= 2;
        }
    }
        
        /* COMPARES DIFFERENT GAP POSITIONS (SHIFTED BY ONE) */
        current_gap_length = true_gap_length;
        if ((count_also_spaced_kmers == 2 || (count_also_spaced_kmers == 1 && current_kmer_length % 2 == 1)) && true_gap_length > 0)
        {
        current_gap_position = true_gap_position + 1;
        lowbit = 1; highbit = 2;
        lowbit <<= ((current_kmer_length - current_gap_position)*2); highbit <<= ((current_kmer_length - current_gap_position)*2);
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));

        /* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH RIGHT BY ONE */
        for(first_half = 0; first_half < 2; first_half++)
        {
        if (current_gap_position < current_kmer_length && current_gap_position > 0 && (count_also_spaced_kmers == 2 ||
        (count_also_spaced_kmers == 1 && ((current_gap_position == current_kmer_length / 2) || current_gap_position == current_kmer_length / 2 + 1))))
        {
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = highbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        }
            
        /* SWITCHES LEFT BY ONE */
        current_gap_position--;
        current_gap_position--;
        compared_kmer = ((current_kmer) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
        lowbit <<= 2; highbit <<= 2;
        }
        }
        
        start = 1;
        end = current_kmer_length;
        
        /* COMPARES UNGAPPED KMER TO ALL SINGLE GAPS IN ALL POSITIONS */
        if (count_also_spaced_kmers != 0 && true_gap_length == 0)
        {
            if(count_also_spaced_kmers == 1)
            {
                start = current_kmer_length / 2;
                end = start + 1 + (current_kmer_length % 2);
            }
        for(current_gap_position = start, current_gap_length = 1; current_gap_position < end; current_gap_position++)
        {
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer << 2) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
            /* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH BEGINNING */
            for(lowbit = 1, highbit = 2, first_half = 0; first_half < 2; first_half++)
            {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
           
        /* END */
        compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
            lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
        }
        }
        }
        
    }
/* Shorter; COMPARES KMER WITH ONE SHORTER */
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
current_kmer_length--;
end = current_kmer_length;
if (current_kmer_length >= shortest_kmer)
{
/* IF NO GAP, INSERTS GAP AT ALL ALLOWED POSITIONS */
if(current_gap_length == 0)
{
    if (count_also_spaced_kmers != 0)
    {
        if(count_also_spaced_kmers == 1)
        {
            start = current_kmer_length / 2;
            end = start + 1 + (current_kmer_length % 2);
        }
        for(current_gap_position = start, current_gap_length = 1; current_gap_position < end; current_gap_position++)
        {
            compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
                if(kmer2_incidence  * kmer_length_difference_cutoff > kmer1_incidence) return (0);
        }
    }
}

current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
    
if (count_also_spaced_kmers != 1) {left = 1; right = 1;}
if (current_gap_position == current_kmer_length / 2 && current_kmer_length % 2 == 0 && count_also_spaced_kmers == 1) {left = 1; right = 0;}
if (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1) left = 1;
    
/* LEFT PART */
if (current_gap_position < current_kmer_length)
{
    if (left == 1 || true_gap_length == 0)
    {
compared_kmer = (current_kmer >> 2);
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
    }
/* RIGHT PART */
if (current_gap_position != 1 && right == 1)
{
if (current_gap_position > 0) current_gap_position--;
    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
current_gap_position = true_gap_position;
/* Shorter with Longer Gap; LONGER GAP */
if(count_also_spaced_kmers != 0 && true_gap_position != 0 && current_gap_length < MAX_GAP_LENGTH)
{
current_gap_length++;
if (current_gap_position < current_kmer_length && left == 1)
{
/* RIGHT BASE GAPPED */
compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    /* printf("DEBUG Localmax: k=%d, gap_pos=%d, gap_len=%d, indices: high=%d, low=%d\n",
           current_kmer_length, current_gap_position, current_gap_length,
           current_kmer_length - current_gap_position,
           current_kmer_length-current_gap_position-1); fflush(stdout); */
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
/* LEFT BASE GAPPED */
if(current_gap_position > 1 && right == 1)
{
current_gap_position--;
compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
current_gap_length--;
current_gap_position++;
}
    /* COMPARES HANGING SINGLE BASE TO UNGAPPED KMER */
    current_gap_position = true_gap_position;
    current_gap_length = true_gap_length;
    if(current_gap_position == 1)
    {
    current_gap_length = 0;
    current_gap_position = 0;
    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
    else if(current_kmer_length-current_gap_position == 0)
    {
    current_gap_length = 0;
    current_gap_position = 0;
    compared_kmer = (current_kmer >> 2);
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
}
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
    
    /* Longer; COMPARES KMER WITH ONE LONGER */
    current_kmer_length++;
    current_kmer_length++;
    if (current_kmer_length < too_long_kmer)
    {
        //if (current_kmer == 0xE34) {  printf("\nDEBUG: TGAnTCA (gap=%d) comparing with longer\n", current_gap_position);}
        compared_kmer = (current_kmer << 2);
        /* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH BEGINNING */
        for(lowbit = 1, highbit = 2, first_half = 0; first_half < 2; first_half++)
        {
            if(count_also_spaced_kmers != 1 || true_gap_length == 0 ||
               (first_half == 1 || (count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2)))
            {
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
                    // if(current_kmer == 0xE34) {
                    // printf("\nDEBUG compare: kmer1 ");
                    // Kmerprint(compared_kmer, current_kmer_length, current_gap_position, current_gap_length);
                    // printf("\tcount=%ld kmer2_count=%ld threshold=%f length=%d gap_pos=%d\n", kmer1_incidence, kmer2_incidence, kmer1_incidence * kmer_length_difference_cutoff, current_kmer_length, current_gap_position);}
                if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
                compared_kmer = lowbit ^ compared_kmer;
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
                
                if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
                compared_kmer = highbit ^ compared_kmer;
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
                
                if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
                compared_kmer = lowbit ^ compared_kmer;
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
                
                if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            }
            /* END */
            compared_kmer = current_kmer;
            lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
            if (true_gap_length == 0) continue;
            current_gap_position++;
            if (current_gap_position > current_kmer_length || (count_also_spaced_kmers == 1 && current_gap_position != (current_kmer_length / 2 + (current_kmer_length % 2)))) break;
        }
        current_gap_position = true_gap_position;
        current_gap_length = true_gap_length;
    

    
/* Longer with Shorter Gap; COMPARES TO ONE SHORTER GAP LENGTH */
if(count_also_spaced_kmers != 0 && true_gap_length >= 1)
{
current_gap_length--;
lowbit = 1; highbit = 2;
lowbit <<= ((current_kmer_length - current_gap_position-1)*2); highbit <<= ((current_kmer_length - current_gap_position-1)*2);
compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
/* NO GAP LEFT */
if (current_gap_length == 0) current_gap_position = 0;
    
    /* LOOP TO ANALYZE ADDED BASE ON EITHER SIDE OF GAP, STARTS WITH RIGHT */
    for(first_half = 0; first_half < 2; first_half++)
    {
        if (current_gap_position < current_kmer_length && (current_gap_position == 0 || Centergap (count_also_spaced_kmers, current_kmer_length, current_gap_position)))
        {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
        }
        
        /* SWITCHES ADDED BASE TO LEFT */
        current_gap_position++;
        if (current_gap_position > current_kmer_length || current_gap_position == 1) break;
        compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    }
current_gap_length++;
}
}
current_gap_position = true_gap_position;
return (1);
}

/* SUBROUTINE THAT GENERATES 128 bit BITMASKS FOR NUCLEOTIDES, KMER STRINGS AND DELETIONS */
void GenerateMask ()
{
    short int counter;
    short int start_position;
    short int position;
    short int current_kmer_length;
    /* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
    for (mask_ULL[1][0] = 3ULL, counter = 0; counter < Nlength; counter++) mask_ULL[1][counter+1] = mask_ULL[1][counter] << 2;
    /* GENERATES mask_ULLS FOR EXTRACTION OF EACH KMER STRING */
    for(current_kmer_length = 2; current_kmer_length <= Nlength; current_kmer_length++)
    {
        for (start_position = 0; start_position < Nlength-current_kmer_length; start_position++)
        {
            for (mask_ULL[current_kmer_length][start_position]=mask_ULL[1][start_position], position = start_position+1; position < current_kmer_length+start_position; position++) {mask_ULL[current_kmer_length][start_position] += mask_ULL[1][position];}
        }
    }
    /* GENERATES HIGH AND LOW mask_ULLS FOR DELETIONS */
    for (lowmask_ULL[0] = mask_ULL[1][0], position = 1; position < Nlength-2; position++) lowmask_ULL[position] = lowmask_ULL[position-1]+mask_ULL[1][position];
    for (highmask_ULL[Nlength-2] = mask_ULL[1][Nlength-2], position = Nlength-3; position > 0; position--) highmask_ULL[position] = highmask_ULL[position+1]+mask_ULL[1][position];
}

// LineBuffer structure to include 5D results
typedef struct {
    uint64_t***** results;  // 5D array [file][kmer_len][gap_pos][gap_len][kmer]
    uint64_t* bitstream;
    uint64_t* temp_bitstream;
    int max_seq_length;
    int file_index;
    char finished;
    pthread_mutex_t mutex;
    pthread_cond_t not_full;
    pthread_cond_t not_empty;
    char** sequences;
    size_t* seq_lengths;
    int write_pos;
    int read_pos;
    int count;
    int shortest_kmer;     // Add these
    int longest_kmer;      // two fields
    FILE* file;
} LineBuffer;


typedef struct PWMConsumerArgs {
    LineBuffer* buffer;
    struct normalized_pwm* signal_pwm;        // Input PWM for finding matches
    struct normalized_pwm* background_pwm;     // Input background PWM for corrections
    struct count_pwm* signal_output_pwm;      // Output PWM for signal counts
    struct count_pwm* background_output_pwm;   // Output PWM for background counts
    double cutoff;
} PWMConsumerArgs;

// Function declarations for externally defined functions
// Function to encode a DNA character to 2-bit representation.
uint64_t encode_base(char base) {
    switch (base) {
        case 'A': return 0x0ULL; // 00
        case 'C': return 0x1ULL; // 01
        case 'G': return 0x2ULL; // 10
        case 'T': return 0x3ULL; // 11
        default: {/*printf("\nError in DNA sequence");*/return UINT64_MAX;} // Return max if invalid character
    }
}

// Function to convert a kmer index back to a DNA string.
void decode_kmer(uint64_t kmer, size_t kmer_length, size_t gap_position, size_t gap_length, char *out) {
    int Ncount = gap_length;
    int kmer_position = kmer_length;
    int i = kmer_position + gap_length - 1;

    for (; i >= 0; i--) {
        if (kmer_position == gap_position && Ncount > 0) {
            out[i] = 'n';
            Ncount--;
        }
        else {
            switch (kmer & 0x3) {
                case 0x0: out[i] = 'A'; break;
                case 0x1: out[i] = 'C'; break;
                case 0x2: out[i] = 'G'; break;
                case 0x3: out[i] = 'T'; break;
            }
            kmer >>= 2;
            kmer_position--;
        }
    }
    out[kmer_length+gap_length] = '\0'; // Null-terminate the string.
}

void print_bitstream_as_dna(const uint64_t *bitstream, size_t block_index, size_t dna_length) {
    printf("Block %zu: ", block_index);
    for (int i = 31; i >= 0; i--) {
        uint64_t mask = 0x3ULL << (2 * i);
        char base = ((32 - i) < dna_length) ?
                    "ACGT"[(bitstream[block_index] & mask) >> (2 * i)] :
                    'a';
        printf("%c", base);
    }
    printf("\n");
}


short int encode_dna_sequence(const char *dna_sequence, uint64_t *bitstream, size_t *bitstream_length) {
    size_t dna_length = strlen(dna_sequence);
    *bitstream_length = (2 * dna_length + 63) / 64;
    memset(bitstream, 0, (*bitstream_length) * sizeof(uint64_t));
    
    // printf("DEBUG: Encoding sequence of length %zu\n", dna_length);
    
    // Initial encoding
    size_t block;
    for (size_t i = 0; i < dna_length; i++) {
        uint64_t encoded_base = encode_base(dna_sequence[i]);
        if (encoded_base == UINT64_MAX) return (0); // DEFECTIVE SEQUENCE
        block = (*bitstream_length) - 1 - i * 2 / 64;
        size_t offset = 62 - (i * 2 % 64);
        bitstream[block] |= (encoded_base << offset);
    }
    bitstream[(*bitstream_length)] = 0;  // add 0 block after so that shifting last block will not bring in garbage seq

    // printf("\nAfter initial encoding:\n");
    // for (size_t i = 0; i < *bitstream_length; i++) {
        // print_bitstream_as_dna(bitstream, i, dna_length);
    // }

    // Handle short sequences
    if (dna_length <= 32) {
        bitstream[0] = bitstream[0] >> (64 - 2 * dna_length);
        // printf("\nAfter short sequence adjustment:\n");
        // print_bitstream_as_dna(bitstream, 0, dna_length);
    }
    
    // Handle longer sequences
    if (dna_length > 32 && dna_length % 32 != 0) {
        size_t left_shift = (dna_length * 2) % 64;
        size_t right_shift = 64 - left_shift;
        
        /*printf("\nBefore shifting blocks:\n");
        for (size_t i = 0; i < *bitstream_length; i++) {
            print_bitstream_as_dna(bitstream, i, dna_length);
            }*/
        
        uint64_t low_bits;
        //Shift blocks maintaining proper bit transitions
        for (size_t block = 0; block < (*bitstream_length); block++) {
            // Get high bits from next block that need to be brought into this block
            uint64_t high_bits = bitstream[block + 1] << left_shift;
            low_bits = (bitstream[block] >> right_shift);
            // Shift current block right and OR in the high bits from next block
            bitstream[block] = low_bits | high_bits;
        }
        
        /*printf("\nAfter shifting blocks:\n");
        for (size_t i = 0; i < *bitstream_length; i++) {
            print_bitstream_as_dna(bitstream, i, dna_length);
        }*/
        
    }
    return(1);
}

short int reverse_complement_bitstream(uint64_t *bitstream, uint64_t *temp_bitstream,
                                     size_t sequence_length) {
    
    size_t bitstream_length = (2 * sequence_length + 63) / 64;
    
    // Copy input to temp_bitstream first
    for(int i = 0; i < bitstream_length; i++) temp_bitstream[i] = bitstream[i];
 
    // printf("\nInput:\n"); for (size_t i = 0; i < bitstream_length; i++) {print_bitstream_as_dna(bitstream, i, sequence_length);}
    
    // Clear original bitstream for reverse complement
    memset(bitstream, 0, bitstream_length * sizeof(uint64_t));
    
    // Process each base pair (2 bits)
    for (size_t i = 0; i < sequence_length; i++) {
        // Find block and position in input (temp_bitstream)
        size_t in_block = (sequence_length - 1 - i) * 2 / 64;
        size_t in_offset = ((sequence_length - 1 - i) * 2) % 64;
        
        // Get 2 bits from temp_bitstream
        uint64_t base = (temp_bitstream[in_block] >> in_offset) & 0x3;
        
        // Complement the base (00->11, 01->10, 10->01, 11->00)
        base = (~base) & 0x3;
        
        // Calculate output position in original bitstream
        size_t out_block = i * 2 / 64;
        size_t out_offset = (i * 2) % 64;
        
        // Write to original bitstream
        bitstream[out_block] |= (base << out_offset);
    }
    
    //printf("\nOutput:\n"); for (size_t i = 0; i < bitstream_length; i++) {print_bitstream_as_dna(bitstream, i, sequence_length);}
    
    return 1;
}

void count_kmers(uint64_t *results, uint64_t *bitstream, size_t sequence_length, size_t kmer_length, size_t gap_position, size_t gap_length) {
    uint64_t main_reg = 0;
    uint64_t source_reg = 0;
    const int shift_amount = 2;
    size_t total_kmers = sequence_length - kmer_length + 1 - gap_length;
    size_t counted_kmers = 0;
    size_t incoming_bits;
    size_t last_source_block_moved = 0;
    uint64_t kmer;
    size_t current_block = 0;
    // char* kmer_str = (char*)malloc(kmer_length + MAX_GAP_LENGTH + 1);
    // kmer_str[kmer_length] = '\0';
    
    // printf("DEBUG: Starting kmer counting, expecting %zu kmers\n", total_kmers);
    
    size_t bitstream_length = (sequence_length * 2 + 63) / 64;

    uint64_t *current_bitstream = bitstream;

    // Load initial blocks
    if (bitstream_length > 0) {
        main_reg = *current_bitstream++;
        bitstream_length--;
        // printf("DEBUG: Loaded block %zu to main: %016lx\n", current_block++, main_reg);
    }
    if (bitstream_length > 0) {
        source_reg = *current_bitstream++;
        bitstream_length--;
        // printf("DEBUG: Loaded block %zu to source: %016lx\n", current_block++, source_reg);
    }

    size_t total_blocks = (sequence_length + 31) / 32;  // How many blocks we should have
    size_t total_number_of_blocks_left_to_move_to_source = total_blocks - 2;
    size_t remaining_bits_in_main = 64; // counts bits left in source and finally main reg

//asm("kmer_count_loop_start:");
    while (counted_kmers < total_kmers) {
        // printf("\nDEBUG: --- Iteration %zu ---\n", counted_kmers);
        // printf("main_reg: %016lx, source_reg: %016lx\n", main_reg, source_reg);
        // printf("remaining_bits_in_main: %zu\n", remaining_bits_in_main);

        // Extract and count kmer
        if (gap_length == 0) {
            kmer = main_reg & ((1ULL << (2 * kmer_length)) - 1);
        } else {
            uint64_t mask = (1ULL << (2 * (kmer_length - gap_position))) - 1;
            kmer = ((main_reg & mask) ^
                   ((main_reg >> (2 * gap_length)) & ~mask)) &
                   ((1ULL << (2 * kmer_length)) - 1);
        }

        // decode_kmer(kmer, kmer_length, gap_position, gap_length, kmer_str);
        // printf("DEBUG: Kmer %zu = %s from block %zu\n", counted_kmers, kmer_str, current_block - 1);
                             
        results[kmer]++;
        counted_kmers++;

        // Shift bits
        uint64_t new_bits = (source_reg & 0x3) << (64 - shift_amount);
        // printf("DEBUG: new_bits: %016lx\n", new_bits);

        #ifdef __x86_64__
        asm volatile (
            "shrdq $2, %1, %0"
            : "+r" (main_reg)
            : "r" (source_reg)
            : "cc"
        );
        #else
        main_reg = (main_reg >> shift_amount) | new_bits;
        #endif

        source_reg >>= shift_amount;
        remaining_bits_in_main -= 2;

        // printf("DEBUG: After shift - main_reg: %016lx, source_reg: %016lx, blocks remaining:%li\n", main_reg, source_reg, total_number_of_blocks_remaining);

        if (last_source_block_moved == 0 && remaining_bits_in_main == 0 && total_number_of_blocks_left_to_move_to_source == 0) {
            remaining_bits_in_main = (sequence_length % 32) * 2;
            last_source_block_moved = 1;
            // printf("DEBUG: last block in source. Remaining blocks, bits %li, %li\n", total_number_of_blocks_left_to_move_to_source, remaining_bits_in_main);
        }
        
        // Critical part: loading next block
        if (last_source_block_moved == 0 && remaining_bits_in_main == 0 && total_number_of_blocks_left_to_move_to_source > 0) {
                source_reg = *current_bitstream++;
                total_number_of_blocks_left_to_move_to_source--;
                remaining_bits_in_main = 64;
                current_block++;
                // printf("DEBUG: new block moved to source. Remaining blocks %li, %li\n", total_number_of_blocks_left_to_move_to_source, remaining_bits_in_main);
            }

    }
//asm("kmer_count_loop_end:");
    
    // free(kmer_str);
    // printf("DEBUG: Counted %zu kmers\n", counted_kmers);
}

uint64_t***** allocate_5d_results(int shortest_kmer, int longest_kmer) {
    uint64_t***** results = malloc(NUM_FILES * sizeof(uint64_t****));
    
    for (int f = 0; f < NUM_FILES; f++) {
        // Allocate for all kmers in range
        results[f] = malloc((longest_kmer + 1) * sizeof(uint64_t***));
        
        for (int k = shortest_kmer; k <= longest_kmer; k++) {
            if (k <= 0) continue;
            
            // Calculate exact size needed for this kmer length
            uint64_t total_kmers = 1ULL << (2 * k);
            
            // Calculate center positions for this k-mer length
            int center1, center2;
            if (k % 2 == 0) {
                center1 = k / 2;
                center2 = -1;  // No second center for even length
            } else {
                center1 = (k - 1) / 2;
                center2 = center1 + 1;
            }
            
            // Allocate array for gap positions
            results[f][k] = malloc((k + 1) * sizeof(uint64_t**));
            for (int g = 0; g <= k; g++) {
                results[f][k][g] = NULL;  // Initialize all to NULL
            }
            
            // Allocate position 0 for ungapped kmers
            results[f][k][0] = malloc(sizeof(uint64_t*));  // Only need index 0
            results[f][k][0][0] = calloc(total_kmers, sizeof(uint64_t));
            if (!results[f][k][0][0]) {
                fprintf(stderr, "Failed to allocate memory for ungapped kmers at k=%d (size=%llu)\n",
                        k, total_kmers);
                exit(1);
            }
            
            // Allocate first center position
            results[f][k][center1] = malloc((MAX_GAP_LENGTH + 1) * sizeof(uint64_t*));
            for (int l = 1; l <= MAX_GAP_LENGTH; l++) {  // Start from 1 for gaps
                results[f][k][center1][l] = calloc(total_kmers, sizeof(uint64_t));
                if (!results[f][k][center1][l]) {
                    fprintf(stderr, "Failed to allocate memory for gapped kmers at k=%d, center=%d, gap=%d (size=%llu)\n",
                            k, center1, l, total_kmers);
                    exit(1);
                }
            }
            
            // Allocate second center position if needed (odd length)
            if (center2 != -1) {
                results[f][k][center2] = malloc((MAX_GAP_LENGTH + 1) * sizeof(uint64_t*));
                for (int l = 1; l <= MAX_GAP_LENGTH; l++) {  // Start from 1 for gaps
                    results[f][k][center2][l] = calloc(total_kmers, sizeof(uint64_t));
                    if (!results[f][k][center2][l]) {
                        fprintf(stderr, "Failed to allocate memory for gapped kmers at k=%d, center=%d, gap=%d (size=%llu)\n",
                                k, center2, l, total_kmers);
                        exit(1);
                    }
                }
            }
        }
    }
    return results;
}

// Matching free function that only frees allocated memory
void free_5d_results(uint64_t***** results, int shortest_kmer, int longest_kmer) {
    if (!results) return;
    
    for (int f = 0; f < NUM_FILES; f++) {
        if (!results[f]) continue;
        
        for (int k = shortest_kmer; k <= longest_kmer; k++) {
            if (k <= 0 || !results[f][k]) continue;
            
            // Free ungapped kmer storage
            if (results[f][k][0]) {
                free(results[f][k][0][0]);
                free(results[f][k][0]);
            }
            
            // Calculate center positions
            int center1, center2;
            if (k % 2 == 0) {
                center1 = k / 2;
                center2 = -1;
            } else {
                center1 = (k - 1) / 2;
                center2 = center1 + 1;
            }
            
            // Free first center position
            if (results[f][k][center1]) {
                for (int l = 1; l <= MAX_GAP_LENGTH; l++) {
                    free(results[f][k][center1][l]);
                }
                free(results[f][k][center1]);
            }
            
            // Free second center position if it exists
            if (center2 != -1 && results[f][k][center2]) {
                for (int l = 1; l <= MAX_GAP_LENGTH; l++) {
                    free(results[f][k][center2][l]);
                }
                free(results[f][k][center2]);
            }
            
            free(results[f][k]);
        }
        free(results[f]);
    }
    free(results);
}

void init_buffer(LineBuffer* buffer, int shortest_kmer, int longest_kmer, int file_index, const char* filename) {
    buffer->file_index = file_index;
    buffer->finished = 0;
    buffer->write_pos = 0;
    buffer->read_pos = 0;
    buffer->count = 0;
    buffer->shortest_kmer = shortest_kmer;
    buffer->longest_kmer = longest_kmer;
    
    pthread_mutex_init(&buffer->mutex, NULL);
    pthread_cond_init(&buffer->not_full, NULL);
    pthread_cond_init(&buffer->not_empty, NULL);
    
    buffer->sequences = malloc(BUFFER_SIZE * sizeof(char*));
    buffer->seq_lengths = malloc(BUFFER_SIZE * sizeof(size_t));
    for (int i = 0; i < BUFFER_SIZE; i++) {
        buffer->sequences[i] = malloc(1024 * sizeof(char));  // Adjust size as needed
    }
    
    buffer->file = fopen(filename, "r");
    if (!buffer->file) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(1);
    }
}

void destroy_buffer(LineBuffer* buffer) {
    for (int i = 0; i < BUFFER_SIZE; i++) {
        free(buffer->sequences[i]);
    }
    free(buffer->sequences);
    free(buffer->seq_lengths);
    
    pthread_mutex_destroy(&buffer->mutex);
    pthread_cond_destroy(&buffer->not_full);
    pthread_cond_destroy(&buffer->not_empty);
    
    fclose(buffer->file);
}

// PRODUCER TO GENERATE LINE BUFFER FOR CONSUMERS
void* producer(void* arg) {
    LineBuffer* buffer = (LineBuffer*)arg;
    char line[1024];
    time_t last_print = time(NULL);
    
    while (fgets(line, sizeof(line), buffer->file)) {
        // Skip empty lines or lines with just whitespace
        size_t len = strlen(line);
        if (len == 0) continue;
        
        // Trim trailing whitespace and newlines
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' || line[len-1] == ' ')) {
            line[--len] = '\0';
        }
        
        // Skip if line is empty after trimming
        if (len == 0) continue;
        
        pthread_mutex_lock(&buffer->mutex);
        
        while (buffer->count == BUFFER_SIZE) {
            time_t current = time(NULL);
            if (current - last_print >= 10) {
                // printf("Producer: Still waiting for buffer space after %ld seconds. Buffer count=%d\n",
                //       current - last_print, buffer->count);
                // fflush(stdout);
                last_print = current;
            }
            pthread_cond_wait(&buffer->not_full, &buffer->mutex);
        }
        
        strcpy(buffer->sequences[buffer->write_pos], line);
        buffer->seq_lengths[buffer->write_pos] = len;
        buffer->write_pos = (buffer->write_pos + 1) % BUFFER_SIZE;
        buffer->count++;
        
        pthread_cond_signal(&buffer->not_empty);
        pthread_mutex_unlock(&buffer->mutex);
    }
    
    pthread_mutex_lock(&buffer->mutex);
    buffer->finished = 1;
    pthread_cond_broadcast(&buffer->not_empty);
    pthread_mutex_unlock(&buffer->mutex);
    
    return NULL;
}

// CONSUMER FOR PWM GENERATION
void* PWM_generator_consumer(void* arg) {
   PWMConsumerArgs* args = (PWMConsumerArgs*)arg;
   LineBuffer* buffer = args->buffer;
   size_t line_count = 1;
   struct match matches;
   struct count_pwm count_pwms[2];
   
   match_init(&matches, args->signal_pwm->width);
   for (int i = 0; i < 2; i++) {
       count_pwm_init(&count_pwms[i], args->signal_pwm->name, args->signal_pwm->width, 0.0);
   }

   // Calculate proper allocation size based on MAX_SEQ_LENGTH
   size_t num_blocks = (MAX_SEQ_LEN + 31) / 32;
   
   // Allocate thread-local bitstreams
   uint64_t* local_bitstream = (uint64_t*)calloc(num_blocks, sizeof(uint64_t));
   uint64_t* local_temp_bitstream = (uint64_t*)calloc(num_blocks, sizeof(uint64_t));
   if (!local_bitstream || !local_temp_bitstream) {
       fprintf(stderr, "Failed to allocate thread-local bitstreams\n");
       free(local_bitstream);
       free(local_temp_bitstream);
       return NULL;
   }

   while (1) {
       pthread_mutex_lock(&buffer->mutex);
       
       while (buffer->count == 0) {
           if (buffer->finished) {
               pthread_mutex_unlock(&buffer->mutex);
               free(local_bitstream);
               free(local_temp_bitstream);
               return NULL;
           }
           pthread_cond_wait(&buffer->not_empty, &buffer->mutex);
           continue;  // Add this to recheck count after waking
       }
       
       char* sequence = buffer->sequences[buffer->read_pos];
       size_t seq_length = buffer->seq_lengths[buffer->read_pos];
       
       // Add validation HERE, before any processing
       if (!sequence || seq_length == 0) {
           printf("\nSkipping empty sequence at position %d", buffer->read_pos);
           pthread_mutex_unlock(&buffer->mutex);
           line_count++;
           continue;
       }
       
       //printf("\nDebug - Processing sequence (len=%zu): '%s'", seq_length, sequence);
       
       buffer->read_pos = (buffer->read_pos + 1) % BUFFER_SIZE;
       buffer->count--;
       
       pthread_cond_signal(&buffer->not_full);
       
       size_t bitstream_length;
       memset(local_bitstream, 0, num_blocks * sizeof(uint64_t));
       
       if (encode_dna_sequence(sequence, local_bitstream, &bitstream_length) == 0) {
           printf("\nError in DNA sequence on line %li, sequence rejected", line_count);
           pthread_mutex_unlock(&buffer->mutex);
           line_count++;
           continue;
       }
       
       pthread_mutex_unlock(&buffer->mutex);
       line_count++;

       if (Findpwmmatch_from_bitstream(args->signal_pwm, args->cutoff,
                                     local_bitstream, seq_length, &matches)) {
           for (int i = 1; i <= matches.position[0]; i++) {
               struct count_pwm* output_pwm = (buffer->file_index == 0) ?
                   args->background_output_pwm : args->signal_output_pwm;
               
               Multinomial_add_to_pwm_from_bitstream(
                   output_pwm,
                   args->signal_pwm,
                   matches.position[i],
                   matches.score[i],
                   args->cutoff,
                   local_bitstream,
                   local_temp_bitstream,
                   seq_length,
                   args->background_pwm);
           }

           reverse_complement_bitstream(local_bitstream, local_temp_bitstream, seq_length);
       
           if (Findpwmmatch_from_bitstream(args->signal_pwm, args->cutoff,
                                         local_bitstream, seq_length, &matches)) {
               for (int i = 1; i <= matches.position[0]; i++) {
                   struct count_pwm* output_pwm = (buffer->file_index == 0) ?
                       args->background_output_pwm : args->signal_output_pwm;
                   
                   Multinomial_add_to_pwm_from_bitstream(
                       output_pwm,
                       args->signal_pwm,
                       matches.position[i],
                       matches.score[i],
                       args->cutoff,
                       local_bitstream,
                       local_temp_bitstream,
                       seq_length,
                       args->background_pwm);
               }
           }
       }
   }
   
   free(local_bitstream);
   free(local_temp_bitstream);
   return NULL;
}

// CONSUMER FOR KMER COUNTING
void* consumer(void* arg) {
    LineBuffer* buffer = (LineBuffer*)arg;
    size_t line_count = 1;
    
    // printf("Consumer: Starting\n");
    while (1) {
        pthread_mutex_lock(&buffer->mutex);
        
        while (buffer->count == 0) {
            if (buffer->finished) {
                pthread_mutex_unlock(&buffer->mutex);
                return NULL;
            }
            pthread_cond_wait(&buffer->not_empty, &buffer->mutex);
        }
        
        char* sequence = buffer->sequences[buffer->read_pos];
        size_t seq_length = buffer->seq_lengths[buffer->read_pos];
        
        buffer->read_pos = (buffer->read_pos + 1) % BUFFER_SIZE;
        buffer->count--;
        
        pthread_cond_signal(&buffer->not_full);

        // Keep mutex locked during sequence processing
        size_t bitstream_length;
        memset(buffer->bitstream, 0, 64 * sizeof(uint64_t));
        
        if (encode_dna_sequence(sequence, buffer->bitstream, &bitstream_length) == 0)
        {
            printf("\nError in DNA sequence on line %li, sequence rejected", line_count);
            pthread_mutex_unlock(&buffer->mutex);
            line_count++;
            continue;
        }
        
        pthread_mutex_unlock(&buffer->mutex);
        line_count++;

        // COUNT KMERS
        for(int strand = 0; strand < 2; strand++)
        {
            
            
            // Process ungapped kmers within valid range
            for (int k = buffer->shortest_kmer; k <= buffer->longest_kmer; k++) {
                if (k <= 0) continue;
                if (buffer->results[buffer->file_index][k][0] &&
                    buffer->results[buffer->file_index][k][0][0]) {
                    count_kmers(buffer->results[buffer->file_index][k][0][0],
                                buffer->bitstream, seq_length, k, 0, 0);
                }
            }
            
            // Process kmers within valid gap position and gap length range
            for (int k = buffer->shortest_kmer; k <= buffer->longest_kmer; k++) {
                if (k <= 0) continue;
                
                // printf("DEBUG: Processing k=%d\n", k);
                
                // Calculate gap positions based on kmer length
                int number_of_gap_positions;
                int gap_positions[2];
                
                if (k % 2 == 0) {
                    number_of_gap_positions = 1;
                    gap_positions[0] = k/2;
                    // printf("DEBUG: Even length kmer, gap position at %d\n", gap_positions[0]);
                } else {
                    number_of_gap_positions = 2;
                    gap_positions[0] = k/2;
                    gap_positions[1] = k/2 + 1;
                    // printf("DEBUG: Odd length kmer, gap positions at %d and %d\n",
                    //        gap_positions[0], gap_positions[1]);
                }
                
                // For each gap position
                for (int i = 0; i < number_of_gap_positions; i++) {
                    int gap_pos = gap_positions[i];
                    // printf("DEBUG: Processing gap position %d\n", gap_pos);
                    
                    // For each gap length (1-10)
                    for (int gap_len = 1; gap_len <= 10; gap_len++) {
                        // Skip if kmer + gap would be longer than sequence
                        if (k + gap_len > seq_length) {
                            // printf("DEBUG: Skipping k=%d gap_len=%d (total %d > seq_length %zu)\n",
                            //        k, gap_len, k + gap_len, seq_length);
                            continue;
                        }
                        
                        // printf("DEBUG: Processing gap length %d\n", gap_len);
                        
                        if (buffer->results[buffer->file_index][k][gap_pos] &&
                            buffer->results[buffer->file_index][k][gap_pos][gap_len]) {
                            
                            // uint64_t array_size = 1ULL << (2 * k);
                            // printf("DEBUG: Result array size for k=%d: %" PRIu64 " entries (%" PRIu64 " bytes)\n",
                            //        k, array_size, array_size * sizeof(uint64_t));
                            
                            count_kmers(buffer->results[buffer->file_index][k][gap_pos][gap_len],
                                        buffer->bitstream, seq_length, k, gap_pos, gap_len);
                            
                            // printf("DEBUG: Finished counting for k=%d gap_pos=%d gap_len=%d\n",
                            //        k, gap_pos, gap_len);
                        }
                    }
                }
            }
            reverse_complement_bitstream(buffer->bitstream, buffer->temp_bitstream, seq_length);
        }
        // printf("Consumer: Finished processing sequence\n\n");
    }
}

// And modify print function to take file parameter:
void print_kmers_single_config(uint64_t***** results, int file, int kmer_len,
                             int shortest_kmer, int longest_kmer,
                             int gap_pos, int gap_len, double threshold,
                             int only_local_max, double length_diff_cutoff) {
    char* kmer_str = malloc(kmer_len + MAX_GAP_LENGTH + 1);
    if (!kmer_str) return;
    

    // if (gap_pos == 0 && gap_len == 0) {printf("UNGAPPED_CHECK k=%d first_kmer_count=%" PRIu64 "\n", kmer_len, results[file][kmer_len][0][0][0]);}
    
    uint64_t total_kmers = 1ULL << (2 * kmer_len);
    
    // printf("\nResults for k=%d gap_pos=%d gap_len=%d:\n",
    //       kmer_len, gap_pos, gap_len);
    
    for (uint64_t kmer = 0; kmer < total_kmers; kmer++) {
        uint64_t signal_count = results[1][kmer_len][gap_pos][gap_len][kmer];
        uint64_t background_count = results[0][kmer_len][gap_pos][gap_len][kmer];

        if (signal_count >= threshold) {
            short int is_localmax = Localmax((long int*****)results,
                                           1, kmer_len,
                                           shortest_kmer, longest_kmer,
                                           gap_pos, gap_len, kmer,
                                           length_diff_cutoff);
            decode_kmer(kmer, kmer_len, gap_pos, gap_len, kmer_str);
            if (!only_local_max || is_localmax == 1) printf("%s\t%i\t%i\t%" PRIu64 "\t%" PRIu64, kmer_str, gap_pos, gap_len, background_count, signal_count);
            if (is_localmax) printf ("\tLocalmax");
            if (!only_local_max || is_localmax == 1) printf("\n");
            }
        }
    free(kmer_str);
}

int main(int argc, char* argv[]) {
    
    
    clock_t start = clock();
    
    short int max_seed_length = 20;
    short int count_kmers = 0;
    uint64_t***** results;
    double threshold;
    double length_diff_cutoff;
    int shortest_kmer;
    int longest_kmer;
    int print_only_localmax;
    struct count_pwm *multinomial_motif;
    short int multinomial;
    double lambda;
    char seed[max_seed_length+1];
    
    struct normalized_pwm qp;
    normalized_pwm_init(&qp, "pwm_for_multinomial", max_seed_length, 0);
    struct count_pwm count_pwms[2];
    
    if (argc != 4 && argc != 5 && argc != 6 && argc != 7) {
        fprintf(stderr, "Usage for kmer counting and local maxima detection: \n\t%s <background_file> <signal_file> <shortest_kmer_length> <longest_kmer_length> [count threshold for printing] [length_diff_cutoff]\n", argv[0]);
        fprintf(stderr, "\tlength_diff_cutoff: Optional parameter for setting the length difference cutoff for local maxima detection (default: 0.35 = longer kmer must have at least 35%% of the counts of the shorter one to be local max)\n");
        fprintf(stderr, "\nUsage for generating a PWM motif using a seed: \n\t%s <background_file> <signal_file> <IUPAC seed> <multinomial> [lambda]\n", argv[0]);
        fprintf(stderr, "\tmultinomial sets how many mismatches are allowed in kmers that contribute to motif, lambda defines non-specific carryover for background correction (not implemented in this version, motif shown is from signal, all matches, uncorrected)\n");
        return 1;
    }


    char* filenames[NUM_FILES] = {argv[1], argv[2]};
    shortest_kmer = atoi(argv[3]);
    //printf("\nkmer length %i", base_kmer_length);
    
    // IF THIRD ARGUMENT IS A STRING, PARSES MOTIF GENERATION ARGUMENTS
    if (shortest_kmer == 0) {
        multinomial = atoi(argv[4]);
        lambda = atof(argv[5]);
        strncpy(seed, argv[3], 20);
        Iupac_to_pwm(&qp, seed);
        //printf("\nMotif generation started with multinomial %i and lambda = %f", multinomial, lambda);
        //fflush(stdout);
        
        // Initialize the output PWMs
        for (int i = 0; i < 2; i++) {
            count_pwm_init(&count_pwms[i], qp.name, qp.width + 2 * FLANK_WIDTH + 1, 0.0);
        }
    }
    else {
        // ELSE PARSES LOCAL MAX ARGUMENTS AND COUNTS KMERS
        count_kmers = 1;
        longest_kmer = atoi(argv[4]);
        threshold = (argc == 6 || argc == 7) ? atof(argv[5]) : 10;             // Use command line value if provided
        length_diff_cutoff = argc == 7 ? atof(argv[6]) : 0.35;  // Use command line value if provided
        
        print_only_localmax = 1;
        // Allocate 5D results array
        results = allocate_5d_results(shortest_kmer, longest_kmer);
        if (!results) {
            fprintf(stderr, "ERROR: Memory allocation failed for results array\n");
            return 1;
        }
    }
    
    // Print parameters
    // printf("DEBUG: Running with parameters:\n");
    // printf("DEBUG: Base kmer length: %d\n", base_kmer_length);
    // printf("DEBUG: Length difference cutoff: %f\n", length_diff_cutoff);
    // printf("DEBUG: Will analyze lengths %d, %d, and %d\n",
    //        base_kmer_length-1, base_kmer_length, base_kmer_length+1);
    
    GenerateMask ();
    
    // Process each file
    for (int file_idx = 0; file_idx < NUM_FILES; file_idx++) {
        LineBuffer buffer;
        buffer.results = results;
        buffer.bitstream = (uint64_t*)calloc(MAX_SEQ_LEN * 2, sizeof(uint64_t));
        buffer.temp_bitstream = (uint64_t*)calloc(MAX_SEQ_LEN * 2, sizeof(uint64_t));
        
        if (!buffer.bitstream || !buffer.temp_bitstream) {
            fprintf(stderr, "ERROR: Memory allocation failed for bitstream\n");
            free_5d_results(results, shortest_kmer, longest_kmer);  // Fixed call
            return 1;
        }

        init_buffer(&buffer, shortest_kmer, longest_kmer, file_idx, filenames[file_idx]);

        pthread_t producer_thread, consumer_thread;
        pthread_create(&producer_thread, NULL, producer, &buffer);
        if (count_kmers == 1) pthread_create(&consumer_thread, NULL, consumer, &buffer);
        else {
            // Create and initialize the PWMConsumerArgs structure
            PWMConsumerArgs pwm_args;
            pwm_args.buffer = &buffer;
            pwm_args.signal_pwm = &qp;
            pwm_args.background_pwm = NULL;  // Or appropriate background PWM if you have one
            pwm_args.signal_output_pwm = &count_pwms[1];
            pwm_args.background_output_pwm = &count_pwms[0];
            pwm_args.cutoff = 0 - (double) multinomial - 0.0001;
            
            pthread_create(&consumer_thread, NULL, PWM_generator_consumer, &pwm_args);
        }

        pthread_join(producer_thread, NULL);
        pthread_join(consumer_thread, NULL);

        destroy_buffer(&buffer);
        free(buffer.bitstream);
    }
    
    if (count_kmers == 1) {
    int file = 0;
    // PRINT KMERS
    // for (int file = 0; file < NUM_FILES; file++) {
    //    printf("\nFile %d results:\n", file);
        printf("kmer\tgap position\tgap length\tbackground\tsignal\tLocalmax\n");
        
        // First print ungapped kmers above threshold for all lengths
        for (int klen = shortest_kmer; klen <= longest_kmer; klen++) {
            if (klen <= 0) continue;
            print_kmers_single_config(results, file, klen, shortest_kmer, longest_kmer,
                                    0, 0, threshold, print_only_localmax, length_diff_cutoff);
        }
        
        // Then print gapped local maxima
        for (int klen = shortest_kmer; klen <= longest_kmer; klen++) {
            if (klen <= 0) continue;
            
            // Calculate center positions
            if (klen % 2 == 0) {
                // Even length: one center position
                int gap_pos = klen/2;
                for (int gap_len = 1; gap_len <= 10; gap_len++) {
                    print_kmers_single_config(results, file, klen, shortest_kmer, longest_kmer,
                                            gap_pos, gap_len, threshold, print_only_localmax, length_diff_cutoff);
                }
            } else {
                // Odd length: two center positions
                for (int gap_pos = klen/2; gap_pos <= klen/2 + 1; gap_pos++) {
                    for (int gap_len = 1; gap_len <= 10; gap_len++) {
                        print_kmers_single_config(results, file, klen, shortest_kmer, longest_kmer,
                                                gap_pos, gap_len, threshold, print_only_localmax, length_diff_cutoff);
                    }
                }
            }
        }
    }
    else {
        int seed_length = strlen(seed);
        printf("\nMotif from all matches to seed in signal file (no background correction)");
        PrintMotif(&count_pwms[1], multinomial, seed_length);  // Print the signal PWM
    }
    
    // Cleanup
    if (count_kmers == 1) free_5d_results(results, shortest_kmer, longest_kmer);
    
    clock_t end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Total execution time: %.2f seconds\n", cpu_time_used);
    
    return 0;
}
