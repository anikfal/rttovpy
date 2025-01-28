def count_lines_between(file_path, start_line, end_line):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Initialize counters and flags
    count = 0
    start_found = False
    for line in lines:
        # Strip whitespace characters from the ends
        stripped_line = line.strip()
        if stripped_line == start_line:
            start_found = True
            continue  # Skip counting the start line itself
        if stripped_line == end_line:
            break  # Stop counting when the end line is found
        if start_found:
            count += 1
    return count