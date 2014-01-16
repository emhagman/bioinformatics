class DNA:

    def __init__(self):
        print('static class, should not see this')

    @staticmethod
    def string_to_list(seq):
        return [s for s in seq]

    @staticmethod
    def rna(dna):
        if isinstance(dna, str):
            dna = dna.upper().replace('T', 'U')
            return DNA.string_to_list(dna)
        else:
            dna = DNA.string_to_list(dna)
            return map(DNA.rna_map, dna)

    @staticmethod
    def complement(seq):
        if isinstance(seq, str):
            seq = DNA.string_to_list(seq)
        return map(DNA.complement_map, seq)

    @staticmethod
    def complement_map(marker):
        marker = marker.upper()
        if 'T' in marker:
            return 'A'
        elif 'A' in marker:
            return 'T'
        elif 'C' in marker:
            return 'G'
        elif 'G' in marker:
            return 'C'

    @staticmethod
    def rna_map(marker):
        marker = marker.upper()
        if 'T' in marker:
            return 'U'
        else:
            return marker
