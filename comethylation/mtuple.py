import re
import itertools



class WithinFragmentComethylationMTuple:
    """A WithinFragmentComethylationMTuple instance stores the within-fragment comethylation counts for a single m-tuple of methylation loci, e.g. a methylation-loci m-tuple.
    
    Attributes:
        chromosome: The chromosome containing the methylation-loci m-tuple.
        chromosome_index: An index to be used when sorting a set of WithinFragmentComethylationMTuple instances by chromosome name.
        m: The "m" in m-tuples, i.e. size of the methylation-loci m-tuples.
        positions: A sorted list of 1-based positions (position_1, position_2, ..., position_m) of length = m-tuple, where position_i is the position of the i-th methylation locus in the methylation-loci m-tuple (with reference to the OT-strand). NB: position_1 < position_2 < ... < position_m by definition.
        methylation_type: The type of methylation event: CG, CHG, CHH, CNN, CG/CHG, CG/CHH, CG/CNN, CHG/CHH, CHG/CNN, CHH/CNN, CG/CHG/CHH, CG/CHG/CNN, CG/CHH/CNN, CHG/CHH/CNN OR CG/CHG/CHH/CNN
        counts: A dictionary storing the counts for each of the 2^m comethylation states combined across strands giving a total of 2^m keys and associated values (counts).
    """
    def __init__(self, chromosome, chromosome_index, m, positions, methylation_type):
        """
        Initiates WithinFragmentComethylationMTuple for a single m-tuple of methylation events with co-ordinates given by arguments (chromosome, positions) and sets all counts to zero.
        """
        if len(positions) != m:
            raise ValueError("__init__ WithinFragmentComethylationMTuple: 'm' must be equal to len(positions)")
        if not all([x in ['CG', 'CHG', 'CHH', 'CNN'] for x in methylation_type.split('/')]):
            raise ValueError("__initi__ WithinFragmentComethylationMTuple: 'methylation_type' must be one or more of 'CG', 'CHG', 'CHH' or 'CNN'. Multiple values must be separated by the '/' character.")

        self.methylation_type = methylation_type
        self.chromosome = chromosome
        self.chromosome_index = chromosome_index
        self.positions = positions
        tmp_k = list((itertools.product(('U', 'M'), repeat = m))) # Step 1 of creating the keys: create all combinations of 'U', 'M' of length m by Cartesian product
        k = [''.join(a) for a in tmp_k] # No need to sort because these are just dictionary keys
        self.counts = dict(zip(k, [0] * (2 ** m))) # Create the dictionary with all counts set to 0
    def display(self): 
        """Display a WithinFragmentComethylationMTuple instance."""
        print 'Methylation type =', self.methylation_type
        print 'Positions =', self.chromosome, ':', self.positions
        print 'Counts =', self.counts
    def m_tuple_id(self):
        print ''.join([self.chromosome, ':', '-'.join([str(a) for a in self.positions])])
    def increment_count(self, comethylation_state, read_1, read_2):
        """Increment the counts attribute based on the comethylation_state that has been extracted from read_1 and read_2. NB: read_2 should be set to None if data is single-end."""
        # Collapse the 9 possible characters in the comethylation_state (., Z, z, X, x, H, h, U, u) to (., M, U, *unchanged*) based on the methylation_type parameter. 
        if 'CG' in self.methylation_type.split('/'):
            comethylation_state = comethylation_state.replace('Z', 'M')
            comethylation_state = comethylation_state.replace('z', 'U')
        if 'CHG' in self.methylation_type.split('/'):
            comethylation_state = comethylation_state.replace('X', 'M')
            comethylation_state = comethylation_state.replace('x', 'U')
        if 'CHH' in self.methylation_type.split('/'):
            comethylation_state = comethylation_state.replace('H', 'M')
            comethylation_state = comethylation_state.replace('h', 'U')
        if 'CNN' in self.methylation_type.split('/'):
            comethylation_state = comethylation_state.replace('U', 'M')
            comethylation_state = comethylation_state.replace('u', 'U')

        # Check whether there are any unexpected, and therefore invalid, characters in comethylation_state
        if not re.search('[^MU]', comethylation_state) is None:
            exit_msg = ''.join([read_1.qname, ' has an invalid comethylation string = ', comethylation_state, '.\nThis should never happen. Please log an issue at www.github.com/PeteHaitch/Comethylation describing the error or email me at peter.hickey@gmail.com.'])
            sys.exit(exit_msg)

        # Single-end
        if read_2 is None and not read_1.is_paired:
            # Check that XG- and XR-tags are compatible directional bisulfite-sequencing protocol. If not then skip the read and report a warning
            if (read_1.opt('XG') == 'CT' and read_1.opt('XR') == 'CT') or (read_1.opt('XG') == 'GA' and read_1.opt('XR') == 'CT'):
                self.counts[comethylation_state] += 1
            else:
                warning_msg = ''.join(['XG-tags or XR-tags for readpair ', read_1.qname, ' are not compatible with the directional bisulfite-sequencing protocol (XG-tags = ', read_1.opt('XG'),', ', read_2.opt('XG'), '; XR-tags = ', read_1.opt('XR'), ', ', read_2.opt('XR'), ')'])
                warnings.warn(warning_msg)
        # Paired-end
        elif read_1.is_paired and read_2.is_paired and read_1.is_read1 and read_2.is_read2:
            # Check that XG- and XR-tags are compatible directional bisulfite-sequencing protocol. If not then skip the read and report a warning
            if (read_1.opt('XG') == 'CT' and read_2.opt('XG') == 'CT' and read_1.opt('XR') == 'CT' and read_2.opt('XR') == 'GA') or (read_1.opt('XG') == 'GA' and read_2.opt('XG') == 'GA' and read_1.opt('XR') == 'CT' and read_2.opt('XR') == 'GA'):
                self.counts[comethylation_state] += 1
            elif not read_1.is_read1 or not read_2.is_read2:
                warning_msg = ''.join(['read_1 or read_2 is incorrectly set for readpair ', read_1.qname])
                warnings.warn(warning_msg)
            else:
                warning_msg = ''.join(['XG-tags or XR-tags for readpair ', read_1.qname, ' are not compatible with the directional bisulfite-sequencing protocol (XG-tags = ', read_1.opt('XG'),', ', read_2.opt('XG'), '; XR-tags = ', read_1.opt('XR'), ', ', read_2.opt('XR'), ')'])
                warnings.warn(warning_msg)             



__all__ = [
    'WithinFragmentComethylationMTuple'
]
