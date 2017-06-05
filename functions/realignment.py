

def smith_waterman(seq1, seq2, score_dict, verbose):

    # The scoring matrix contains an extra row and column for the gap (-), hence
    # the +1 here.
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    # Initialize the scoring matrix.
    score_matrix, start_pos, score = create_score_matrix(rows, cols, seq1, seq2, score_dict)

    # Traceback. Find the optimal path through the scoring matrix. This path
    # corresponds to the optimal local sequence alignment.
    seq1_aligned, seq2_aligned = traceback(score_matrix, start_pos, seq1, seq2)
    assert len(seq1_aligned) == len(seq2_aligned), 'aligned strings are not the same size'

    if verbose:
        # Pretty print the results. The printing follows the format of BLAST results
        # as closely as possible.
        alignment_str, idents, gaps, mismatches = alignment_string(seq1_aligned, seq2_aligned)
        alength = len(seq1_aligned)
        print
        print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,
              alength, idents / alength, gaps, alength, gaps / alength))
        print

        seq1_slice = seq1_aligned
        print('Query  {0:<4}  {1}  {2:<4}'.format(1, seq1_slice, len(seq1_slice)))
        print('             {0}'.format(alignment_str))
        seq2_slice = seq2_aligned
        print('Sbjct  {0:<4}  {1}  {2:<4}'.format(1, seq2_slice, len(seq2_slice)))
        print
        print (start_pos[0] - len(seq2))
        print '#Max Score:', score
    return (start_pos[0] - len(seq2)), score


def create_score_matrix(rows, cols, seq1, seq2, score_dict):
    '''Create a matrix of scores representing trial alignments of the two sequences.

    Sequence alignment can be treated as a graph search problem. This function
    creates a graph (2D matrix) of scores, which are based on trial alignments
    of different base pairs. The path with the highest cummulative score is the
    best alignment.
    '''
    score_matrix = [[0 for col in range(cols)] for row in range(rows)]
    # Fill the scoring matrix.
    max_score = 0
    max_pos   = None    # The row and columbn of the highest score in matrix.
    for i in range(1, rows):
        for j in range(1, cols):
            score = calc_score(score_matrix, i, j, seq1, seq2, score_dict)
            if score > max_score:
                max_score = score
                max_pos   = (i, j)

            score_matrix[i][j] = score

    assert max_pos is not None, 'the x, y position with the highest score was not found'

    return score_matrix, max_pos, max_score

def calc_score(matrix, x, y, seq1, seq2, score_dict):
    '''Calculate score for a given x, y position in the scoring matrix.

    The score is based on the up, left, and upper-left neighbors.
    '''
    match       = score_dict['match']
    mismatch    = score_dict['mismatch']
    gap         = score_dict['gap']

    similarity = match if seq1[x - 1] == seq2[y - 1] else mismatch

    diag_score = matrix[x - 1][y - 1] + similarity
    up_score   = matrix[x - 1][y] + gap
    left_score = matrix[x][y - 1] + gap

    return max(0, diag_score, up_score, left_score)



def traceback(score_matrix, start_pos, seq1, seq2):
    '''Find the optimal path through the matrix.

    This function traces a path from the bottom-right to the top-left corner of
    the scoring matrix. Each move corresponds to a match, mismatch, or gap in one
    or both of the sequences being aligned. Moves are determined by the score of
    three adjacent squares: the upper square, the left square, and the diagonal
    upper-left square.

    WHAT EACH MOVE REPRESENTS
        diagonal: match/mismatch
        up:       gap in sequence 1
        left:     gap in sequence 2
    '''

    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y         = start_pos
    move         = next_move(score_matrix, x, y)
    while move != END:
        if move == DIAG:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append(seq2[y - 1])
            x -= 1
            y -= 1
        elif move == UP:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append('-')
            x -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[y - 1])
            y -= 1

        move = next_move(score_matrix, x, y)

    aligned_seq1.append(seq1[x - 1])
    aligned_seq2.append(seq2[y - 1])

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


def next_move(score_matrix, x, y):
    diag = score_matrix[x - 1][y - 1]
    up   = score_matrix[x - 1][y]
    left = score_matrix[x][y - 1]
    if diag >= up and diag >= left:     # Tie goes to the DIAG move.
        return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
    elif up > diag and up >= left:      # Tie goes to UP move.
        return 2 if up != 0 else 0      # UP move or end.
    elif left > diag and left > up:
        return 3 if left != 0 else 0    # LEFT move or end.
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')


def alignment_string(aligned_seq1, aligned_seq2):
    '''Construct a special string showing identities, gaps, and mismatches.

    This string is printed between the two aligned sequences and shows the
    identities (|), gaps (-), and mismatches (:). As the string is constructed,
    it also counts number of identities, gaps, and mismatches and returns the
    counts along with the alignment string.

    AAGGATGCCTCAAATCGATCT-TTTTCTTGG-
    ::||::::::||:|::::::: |:  :||:|   <-- alignment string
    CTGGTACTTGCAGAGAAGGGGGTA--ATTTGG
    '''
    # Build the string as a list of characters to avoid costly string
    # concatenation.
    idents, gaps, mismatches = 0, 0, 0
    alignment_string = []
    for base1, base2 in zip(aligned_seq1, aligned_seq2):
        if base1 == base2:
            alignment_string.append('|')
            idents += 1
        elif '-' in (base1, base2):
            alignment_string.append(' ')
            gaps += 1
        else:
            alignment_string.append(':')
            mismatches += 1

    return ''.join(alignment_string), idents, gaps, mismatches

def reverse_strand(in_str):
    return_string = ''
    for char in in_str:
        if char == 'a' or char == 'A':
            return_string = 'T' + return_string
        elif char == 't' or char == 'T':
            return_string = 'A' + return_string
        elif char == 'c' or char == 'C':
            return_string = 'G' + return_string
        elif char == 'g' or char == 'G':
            return_string = 'C' + return_string
    return return_string


def expansion_aware_realign(sample, pre, post, motif, score_dict, verbose = False):
    read_len = len(pre)
    max_score = 0
    max_nCopy = 0
    max_pos = 0
    for nCopy in range((int(read_len / len(motif)) + 1) + 1):
        if verbose:
            print '####### nCopy =', nCopy, ' ##########'
        var_realign_string = pre + motif * nCopy + post
        pos, score = smith_waterman(var_realign_string, sample, score_dict, verbose)
        if verbose:
            print 'Score = ', score
            print 'Best Alignment Position = ', pos
        if score > max_score:
            max_score = score
            max_nCopy = nCopy
            max_pos = pos
        if score == read_len * score_dict['match']:
            break
    if verbose:
        print '>>>>>>>'
        print '> Best nCopy = ', max_nCopy
        print '> Score (for best nCopy) = ', max_score
        print '> Best Alignment Position = ', max_pos
    return max_nCopy, max_pos, max_score

# margin: amount of slip we allow between alignment position and STR start and end
def classify_realigned_read(sample, motif, start_pos, nCopy, score, score_dict, read_len, margin, verbose):
    end_pos = start_pos + len(sample) - 1

    start_str = read_len
    end_str = read_len + nCopy * len(motif)

    if verbose:
        print 'STR region: ', [start_str, end_str]
        print 'Realigned Sample region: ', [start_pos, end_pos]


    if start_pos >= start_str - margin and start_pos <= end_str + margin:
        start_in_str = True
    else:
        start_in_str = False

    if end_pos >= start_str - margin and end_pos <= end_str + margin:
        end_in_str = True
    else:
        end_in_str = False

    score_threshold = int(0.9 * read_len * score_dict['match']) 

    if score < score_threshold:
        if verbose:
            print '>> Non-spanning (low score)'
            print '>> Score:', score
        return 'NoSpan'
    else:
        if nCopy == 0:
            if verbose:
                print '>> Non-spanning'
            return 'NoSpan'
        else:
            if start_in_str and end_in_str:
                if verbose:
                    print '>> IRR'
                return 'IRR'
            elif start_in_str and ~end_in_str:
                if verbose:
                    print '>> Post-Flank'
                return 'PostFlank'
            elif ~start_in_str and end_in_str:
                if verbose:
                    print '>> Pre-Flank'
                return 'PreFlank'
            elif ~start_in_str and ~end_in_str:
                if start_pos < start_str and end_pos > end_str:
                    if verbose:
                        print '>> Enclosing'
                    return 'Enclosing'
                else:
                    if verbose:
                        print 'UNKNOWN TYPE'
                    return 'unknown'
