import numpy as np

def smith_waterman_all_alignments(seq1, seq2, match_score=2, mismatch_score=-1, gap_penalty=-2):

    # Initialize scoring matrix
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = np.zeros((rows, cols), dtype=int)

    # Store max score and its position
    max_score = 0
    max_positions = []

    # Fill the scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):

            match = matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            # Calculate scores
            gap_in_seq1 = matrix[i - 1][j] + gap_penalty  # Gap penalty for inserting a gap in seq1
            gap_in_seq2 = matrix[i][j - 1] + gap_penalty  # Gap penalty for inserting a gap in seq2

            # Update matrix with the best score (taking max of all three conditions)
            matrix[i][j] = max(0, match, gap_in_seq1, gap_in_seq2)

            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_positions = [(i, j)]
            elif matrix[i][j] == max_score:
                max_positions.append((i, j))

    print("Scoring Matrix:")
    print(matrix)
    print(f"\nMax Score: {max_score} at positions {[(i + 1, j + 1) for i, j in max_positions]} in matrix\n")

    # Traceback for obtaining best local alignment
    def traceback(start_pos):
        aligned1, aligned2 = [], []
        i, j = start_pos
        while matrix[i][j] != 0:
            score = matrix[i][j]
            diag = matrix[i - 1][j - 1]
            up = matrix[i - 1][j]
            left = matrix[i][j - 1]

            if score == diag + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score):
                aligned1.append(seq1[i - 1])
                aligned2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif score == up + gap_penalty:
                aligned1.append(seq1[i - 1])
                aligned2.append('-')
                i -= 1
            elif score == left + gap_penalty:
                aligned1.append('-')
                aligned2.append(seq2[j - 1])
                j -= 1
            else:
                break

        # Reverse the sequences since traceback goes backwards
        return ''.join(reversed(aligned1)), ''.join(reversed(aligned2))

    # Run traceback from all max positions
    alignments = [traceback(pos) for pos in max_positions]

    # Choose the best alignment (by length)
    best_alignment = max(alignments, key=lambda x: len(x[0]))

    return best_alignment

def get_score_input(prompt, default):
    try:
        user_input = input(f"{prompt} [default {default}]: ")
        if user_input.strip() == "":
            print(f"No input given. Using default {prompt.lower()}: {default}")
            return default
        return int(user_input)
    except ValueError:
        print(f"Invalid input. Using default {prompt.lower()}: {default}")
        return default


def main():
    print(f"# Example sequence (input can be given either in lower case or upper case letters)\n# seq1_str = {"GCGTTACA"}\n# seq2_str = {"GCTTGCA"}\n")
    print("Note: If early or weaker matches are not appearing in the alignment, consider using a higher match score or lower penalties.\n")

    try:
        seq1 = input("Enter Sequence 1: ").upper()
        seq2 = input("Enter Sequence 2: ").upper()

        # Check for empty inputs
        if not seq1 or not seq2:
            raise ValueError("Both sequences must be provided and non-empty.")

        valid_chars = set("ATGCU") # standard accepted nucleotides
        if not (set(seq1).issubset(valid_chars) and set(seq2).issubset(valid_chars)):
            raise ValueError("Invalid characters detected. Allowed characters are A, T, G, C, U only.")

        # Get scoring values with fallback + messages
        match_score = get_score_input("Enter match score", 2)
        mismatch_score = get_score_input("Enter mismatch score", -1)
        gap_penalty = get_score_input("Enter gap penalty", -2)

        # Run alignment
        a1, a2 = smith_waterman_all_alignments(seq1, seq2)
        print("Best Aligned A:", a1)
        print("Best Aligned B:", a2)

    except Exception as e:
        print("Error:", e)

if __name__ == "__main__":
    main()

