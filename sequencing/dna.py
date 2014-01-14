class DNA:

    def __init__(self):
        print('static class, should not see this')

    @staticmethod
    def string_to_list(seq):
        return [s for s in seq]

    @staticmethod
    def mrna(dna):
        if 'string' in type(dna):
            dna = dna.lower().replace('a', 'u')
            return DNA.string_to_list(dna)
        else:
            dna = DNA.string_to_list(dna)
            return map(DNA.mrna_map, dna)

    @staticmethod
    def complement(seq):
        if 'string' in type(seq):
            seq = DNA.string_to_list(seq)
        return map(DNA.complement_map, seq)

    @staticmethod
    def complement_map(marker):
        marker = marker.lower()
        if 't' in marker:
            return 'a'
        elif 'a' in marker:
            return 't'
        elif 'c' in marker:
            return 'g'
        elif 'g' in marker:
            return 'c'

    @staticmethod
    def mrna_map(marker):
        marker = marker.lower()
        if 'a' in marker:
            return 'u'
        else:
            return marker
