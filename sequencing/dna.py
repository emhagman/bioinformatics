class DNA:

    def __init__(self):
        print('static class, should not see this')

    @staticmethod
    def complement(seq):
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