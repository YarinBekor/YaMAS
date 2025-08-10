# This function using certain flags to find the unique speices within the sample
# It won't count any species that don't have files on the NCBI database
def search_species(filename, output_directory, dictionary_path):
    with open(filename, 'r') as file:
        lines = file.readlines()

    comparison_filename = dictionary_path
    with open(comparison_filename, 'r') as comparison_file:
        comparison_lines = comparison_file.readlines()

    comparison_items = set()
    for line in comparison_lines:
        comparison_items.add(line.strip())

    items = set()
    for line in lines:
        start_index = 0
        while True:
            start_index = line.find('|s__', start_index)
            if start_index == -1:
                break
            end_index = line.find('2|', start_index)
            if end_index == -1:
                break
            item = line[start_index + 3: end_index].strip()
            if '|t__' in item:
                item = item.split('|t__')[0]
            if item.startswith('_'):
                item = item[1:]
            if item in comparison_items:
                items.add(item)
            start_index = end_index + 2

    if items:
        max_length = max(len(item) for item in items)
        formatted_items = [item.ljust(max_length) for item in items]

        output_file = f"{output_directory}/species_list.txt"
        with open(output_file, 'w') as output_file:
            for item in formatted_items:
                output_file.write(item + "\n")

    return list(items)

