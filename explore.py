import h5py
import re

def extract_coordinates(h5_file_path, output_file_path):
    """
    Extract genomic coordinates from H5 dataset names and save them to a text file.
    """
    coordinates = []

    # Open the H5 file and extract dataset names
    with h5py.File(h5_file_path, 'r') as f:
        # Get all dataset names
        def collect_names(name, obj):
            if isinstance(obj, h5py.Dataset):
                coordinates.append(name)
        f.visititems(collect_names)

    # Process each coordinate and format for annotation file
    formatted_coords = []
    for coord in coordinates:
        # Parse the coordinates using regex
        match = re.match(r'(chr\d+)_(\d+)_(\d+)_predictions', coord)
        if match:
            chrom, start, end = match.groups()
            # Format: chromosome    start    end    region_name
            formatted_line = f"{chrom}\t{start}\t{end}\t{coord}"
            formatted_coords.append(formatted_line)

    # Write to output file
    with open(output_file_path, 'w') as f:
        # Write header
        f.write("chromosome\tstart\tend\tregion_name\n")
        # Write coordinates
        for line in formatted_coords:
            f.write(line + '\n')

    print(f"Coordinates have been extracted to: {output_file_path}")
    print(f"Total regions extracted: {len(formatted_coords)}")

# Usage example:
extract_coordinates('/Users/gabrielduarte/Documents/GitHub/scPrediXcan/enformer_preds/Step2_Personalized_epigenomics.h5', 'gene_anno2.txt')
