# Written by Jesse Gomer (jessegomer@gmail.com)
# for the University of California Conservation Consortium's CALeDNA Program

# example usage python reformat_summary_for_r.py nochim_mergedCO1.txt  r_sum.txt 95
import argparse
import shutil


def summarize_taxonomy(full_taxonomy, confidences):
    full_taxonomy = full_taxonomy.rstrip(';')
    confidences = confidences.rstrip(';')
    taxonomy_components = [f.split(":") for f in full_taxonomy.split(";")]
    confidence_components = [c.split(":") for c in confidences.split(";")]
    summary = ";".join([":".join((t[0], t[1], confidence_components[i][1])) for i, t in enumerate(taxonomy_components)])
    return summary


def truncate_taxonomy(full_taxonomy, confidences, cutoff):
    # the taxonomy and confidences may have an extra semicolon at the end
    full_taxonomy = full_taxonomy.rstrip(';')
    confidences = confidences.rstrip(';')
    taxonomy = dict([level.split(':', 1) for level in full_taxonomy.split(';')])
    truncated_taxonomy = {}
    for level_info in confidences.split(';'):
        level_name, confidence_value = level_info.split(':')
        if float(confidence_value) >= cutoff:
            truncated_taxonomy[level_name] = taxonomy[level_name]
    return truncated_taxonomy


def reformat_summary(summary_file_name, output_file_name, cutoff):
    output_levels = ["superkingdom","phylum", "class", "order", "family", "genus", "species"]

    summary = open(summary_file_name).readlines()
    previous_header = summary[0].strip().split('\t')
    taxonomy_index = previous_header.index('taxonomy')
    confidence_index = previous_header.index('taxonomy_confidence')
    header = previous_header[:taxonomy_index] + ['sum.taxonomy', 'taxonomy_confidence']
    output = open(output_file_name + '.tmp', 'w')
    output.write('\t'.join(header) + '\n')

    for line in summary[1:]:
        fields = line.strip('\n').split('\t')
        # a colon in the taxonomy means that something was found
        if ':' in fields[taxonomy_index]:
            taxonomy = truncate_taxonomy(fields[taxonomy_index], fields[confidence_index], cutoff)
            output_taxonomy = [taxonomy.get(level, '') for level in output_levels]
            taxonomy_summary = summarize_taxonomy(fields[taxonomy_index], fields[confidence_index])
        else:
            output_taxonomy = ''
            taxonomy_summary = ''
        fields_to_write = fields[:taxonomy_index] + [';'.join(output_taxonomy)] + [taxonomy_summary]
        output.write('\t'.join(fields_to_write) + '\n')

    output.close()
    shutil.move(output_file_name + '.tmp', output_file_name)


parser = argparse.ArgumentParser(description='Reformats a summary table for use in R code')
parser.add_argument('summary_file', type=str,  help='Summary file')
parser.add_argument('output_file', type=str, help='File where output will be written')
parser.add_argument('cutoff', type=float, help='Confidence percent cutoff to include [0-100]')

if __name__ == '__main__':
    args = parser.parse_args()
    summary_file = args.summary_file
    cutoff = args.cutoff
    output_file = args.output_file
    reformat_summary(summary_file, output_file, cutoff)
