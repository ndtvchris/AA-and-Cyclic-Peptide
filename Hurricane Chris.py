class cyclone:
    aminos = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
              'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

    def __init__(self):
        self.peptide = 'HCVCPVDILN'
        self.fragments = {}
        self.pepLen = len(self.peptide)
        self.fragments[' '] = 0
        self.aminoList = []
        self.aminoList.append(' ')

    def chunkyFilas(self):
        window = 1
        holdList = []
        while window <= self.pepLen:
            for i in range(self.pepLen):
                for e in range(i, i + window):
                    try:
                        holdList.append(self.peptide[e])
                    except:
                        holdList.append(self.peptide[(e - self.pepLen)])
                self.fragments[''.join(holdList)] = 0
                self.aminoList.append(''.join(holdList))
                if window == self.pepLen:
                    break
                holdList.clear()
            window += 1

    def bigMasses(self):

        for key in self.fragments:
            for amino in key:
                if amino in cyclone.aminos:
                    self.fragments[key] += cyclone.aminos[amino]

    def lineUp(self):
        pList = []
        for key in self.aminoList:
            pList.append(self.fragments[key])
        pList = sorted(pList)
        pList = [str(i) for i in pList]
        print(' '.join(pList))


c = cyclone()
c.chunkyFilas()
c.bigMasses()
c.lineUp()
