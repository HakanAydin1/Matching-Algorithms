import time
import random
import string
import matplotlib.pyplot as plt


# Brute-Force Algorithm
def brute_force(pattern, text):
    m = len(pattern)
    n = len(text)
    occurrences = []

    for i in range(n - m + 1):
        j = 0
        while j < m and text[i + j] == pattern[j]:
            j += 1

        if j == m:
            occurrences.append(i)

    return occurrences

# Sunday Algorithm
def sunday(pattern, text):
    m = len(pattern)
    n = len(text)
    occurrences = []

    def calculate_shifts():
        shifts = {ch: m + 1 for ch in text}
        for i in range(m):
            shifts[pattern[i]] = m - i
        return shifts

    shifts = calculate_shifts()

    i = 0
    while i <= n - m:
        j = 0
        while j < m and text[i + j] == pattern[j]:
            j += 1

        if j == m:
            occurrences.append(i)

        if i + m < n:
            i += shifts[text[i + m]]
        else:
            break

    return occurrences

# Knuth-Morris-Pratt (KMP) Algorithm
def kmp(pattern, text):
    m = len(pattern)
    n = len(text)
    occurrences = []

    def calculate_lps():
        lps = [0] * m
        i = 1
        length = 0

        while i < m:
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1

        return lps

    lps = calculate_lps()

    i = 0
    j = 0
    while i < n:
        if pattern[j] == text[i]:
            i += 1
            j += 1

            if j == m:
                occurrences.append(i - j)
                j = lps[j - 1]
        else:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1

    return occurrences

# Finite State Machine (FSM) Algorithm
def fsm(text, pattern):
    def build_transition_table(pattern):
        transitions = [0]
        state = 0

        for i in range(1, len(pattern)):
            while state > 0 and pattern[state] != pattern[i]:
                state = transitions[state - 1]
            if pattern[state] == pattern[i]:
                state += 1
            else:
                state = 0
            transitions.append(state)

        return transitions

    transitions = build_transition_table(pattern)
    matches = []

    state = 0
    for i, char in enumerate(text):
        while state > 0 and pattern[state] != char:
            state = transitions[state - 1]
        if pattern[state] == char:
            state += 1
            if state == len(pattern):
                matches.append(i - len(pattern) + 1)
                state = transitions[state - 1]

    return matches


# Rabin-Karp Algorithm
def rabin_karp(pattern, text):
    m = len(pattern)
    n = len(text)
    occurrences = []

    def calculate_hash(string):
        h = 0
        for char in string:
            h = (256 * h + ord(char)) % 101
        return h

    def rehash(prev_hash, prev_char, next_char):
        return (prev_hash * 256 - ord(prev_char) * 256 ** m + ord(next_char)) % 101

    pattern_hash = calculate_hash(pattern)
    text_hash = calculate_hash(text[:m])

    for i in range(n - m + 1):
        if pattern_hash == text_hash:
            if pattern == text[i:i + m]:
                occurrences.append(i)

        if i < n - m:
            text_hash = rehash(text_hash, text[i], text[i + m])

    return occurrences

# Gusfield Z Algorithm
def gusfield_z(pattern, text):
    def calculate_z_array(concatenated):
        n = len(concatenated)
        z = [0] * n
        left = right = 0

        for i in range(1, n):
            if i > right:
                left = right = i
                while right < n and concatenated[right - left] == concatenated[right]:
                    right += 1
                z[i] = right - left
                right -= 1
            else:
                k = i - left
                if z[k] < right - i + 1:
                    z[i] = z[k]
                else:
                    left = i
                    while right < n and concatenated[right - left] == concatenated[right]:
                        right += 1
                    z[i] = right - left
                    right -= 1

        return z

    concatenated = pattern + "$" + text
    z_array = calculate_z_array(concatenated)
    m = len(pattern)
    n = len(text)
    occurrences = []

    for i in range(m + 1, m + n + 1):
        if z_array[i] == m:
            occurrences.append(i - m - 1)

    return occurrences

# Load book chapters as text
def load_book_chapters():
    book_chapters = []
    for i in range(1, 6):  # Assuming there are 5 chapters
        chapter_file = f"book/chapter{i}.txt"  # Replace with actual chapter file names
        with open(chapter_file, "r") as file:
            chapter_text = file.read()
            book_chapters.append(chapter_text)
    return book_chapters

def generate_random_pattern(length):
    pattern = ''.join(random.choice(string.ascii_lowercase) for _ in range(length))
    return pattern

# Generate patterns of different lengths
def generate_patterns(pattern_lengths):
    patterns = []
    for length in pattern_lengths:
        pattern = generate_random_pattern(length)  # Replace with your desired pattern generation logic
        patterns.append(pattern)
    return patterns

# Measure the running times of pattern matching algorithms
def measure_running_times(algorithms, chapters, patterns):
    running_times = {}

    for algorithm in algorithms:
        times = []
        for pattern in patterns:
            start_time = time.time()
            for chapter in chapters:
                algorithm(chapter, pattern)
            end_time = time.time()
            elapsed_time = end_time - start_time
            times.append(elapsed_time)

        running_times[algorithm.__name__] = times

    return running_times


# Plot the running time comparison graph
def plot_running_times(text_lengths, running_times):
    algorithms = list(running_times.keys())

    for algorithm in algorithms:
        times = running_times[algorithm]
        plt.plot(text_lengths[:len(times)], times, label=algorithm)

    plt.xlabel('Text Length')
    plt.ylabel('Running Time (seconds)')
    plt.legend()
    plt.show()

# Task 1: Compare running times of pattern matching algorithms using small and large patterns
def compare_running_times():
    book_chapters = load_book_chapters()

    small_pattern = generate_random_pattern(10)  # Replace with your desired small pattern length
    large_pattern = generate_random_pattern(500)  # Replace with your desired large pattern length

    patterns = [small_pattern, large_pattern]

    algorithms = [brute_force, sunday, kmp, fsm, rabin_karp, gusfield_z]

    running_times = measure_running_times(algorithms, book_chapters, patterns)

    text_lengths = [len(chapter) for chapter in book_chapters]

    plot_running_times(text_lengths, running_times)

# Task 2: Prove empirical results for specific scenarios
def prove_empirical_results():
    T = generate_random_pattern(102400)  # Replace with your desired text length
    P1 = T[:50000]  # Replace with your desired pattern lengths
    P2 = T[:100000]
    P3 = T[:150000]

    algorithms = [sunday, kmp, rabin_karp]

    running_times = measure_running_times(algorithms, [T], [P1, P2, P3])

    pattern_lengths = [len(pattern) for pattern in [P1, P2, P3]]

    plot_running_times(pattern_lengths, running_times)

# Run the tasks
compare_running_times()
prove_empirical_results()