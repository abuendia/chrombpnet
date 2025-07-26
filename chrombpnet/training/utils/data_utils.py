import numpy as np
import pandas as pd
import pyBigWig
import pyfaidx
import pysam
from chrombpnet.training.utils import one_hot
from tqdm import tqdm


class PersonalizedGenome:
    """
    A class to handle personalized genomes by applying variants from VCF to reference genome.
    """
    def __init__(self, reference_genome, vcf_file=None, sample_id=None):
        """
        Initialize personalized genome.
        
        Args:
            reference_genome: pyfaidx.Fasta object of reference genome
            vcf_file: Path to VCF/BCF file containing variants
            sample_id: Sample ID to extract genotypes from VCF
        """
        self.reference_genome = reference_genome
        self.vcf_file = vcf_file
        self.sample_id = sample_id
        self.variants_cache = {}
        
        if vcf_file and sample_id:
            self._load_variants()
    
    def _load_variants(self):
        """Load variants from VCF file for the specified sample."""
        print(f"Loading variants from {self.vcf_file} for sample {self.sample_id}")
        
        with pysam.VariantFile(self.vcf_file, "rb") as vcf:
            if self.sample_id not in vcf.header.samples:
                raise ValueError(f"Sample {self.sample_id} not found in VCF file")
            
            for record in vcf.fetch():
                # Get genotype for the sample
                sample_gt = record.samples[self.sample_id]['GT']
                if sample_gt is None or sample_gt == (None, None):
                    continue
                
                # Fix: Check if ANY allele is alternate (not just the first)
                has_alt = any(gt is not None and gt > 0 for gt in sample_gt)
                if not has_alt:
                    continue
                    
                chrom = record.chrom
                pos = record.pos - 1  # Convert to 0-based
                ref_allele = record.ref
                alt_alleles = record.alts
                
                # Get all alternate alleles for this sample
                alt_indices = [gt for gt in sample_gt if gt is not None and gt > 0]
                
                # For heterozygous variants, we need to handle both alleles
                if len(alt_indices) == 1 and len(set(sample_gt)) == 2:
                    # Heterozygous: 0/1 or 1/0
                    alt_allele = alt_alleles[alt_indices[0] - 1]
                    is_heterozygous = True
                elif len(alt_indices) == 2:
                    # Homozygous alternate: 1/1
                    alt_allele = alt_alleles[alt_indices[0] - 1]  # Use first alt allele
                    is_heterozygous = False
                else:
                    continue
                
                # Store variant information
                if chrom not in self.variants_cache:
                    self.variants_cache[chrom] = []
                self.variants_cache[chrom].append({
                    'pos': pos,
                    'ref': ref_allele,
                    'alt': alt_allele,
                    'is_heterozygous': is_heterozygous,
                    'genotype': sample_gt
                })
        
        # Sort variants by position for each chromosome
        for chrom in self.variants_cache:
            self.variants_cache[chrom].sort(key=lambda x: x['pos'])
        
        print(f"Loaded {sum(len(vars) for vars in self.variants_cache.values())} variants")
    
    def get_sequence(self, chrom, start, end):
        """
        Get personalized sequence for a genomic region.
        
        Args:
            chrom: Chromosome name
            start: Start position (0-based)
            end: End position (0-based)
        
        Returns:
            Personalized DNA sequence string
        """
        # Get reference sequence
        ref_seq = str(self.reference_genome[chrom][start:end])
        
        if not self.vcf_file or chrom not in self.variants_cache:
            return ref_seq
        
        # Apply variants in the region
        personalized_seq = list(ref_seq)
        
        for variant in self.variants_cache[chrom]:
            var_pos = variant['pos']
            if start <= var_pos < end:
                # Adjust position relative to our region
                rel_pos = var_pos - start
                
                # Apply the variant
                ref_allele = variant['ref']
                alt_allele = variant['alt']
                
                # Check if reference allele matches
                if rel_pos + len(ref_allele) <= len(personalized_seq):
                    region_ref = ''.join(personalized_seq[rel_pos:rel_pos + len(ref_allele)])
                    if region_ref == ref_allele:
                        # Replace reference with alternate
                        if len(alt_allele) == len(ref_allele):
                            # Simple substitution
                            for i, base in enumerate(alt_allele):
                                personalized_seq[rel_pos + i] = base
                        elif len(alt_allele) > len(ref_allele):
                            # Insertion
                            for i, base in enumerate(alt_allele):
                                if i < len(ref_allele):
                                    personalized_seq[rel_pos + i] = base
                                else:
                                    personalized_seq.insert(rel_pos + len(ref_allele), base)
                        else:
                            # Deletion
                            for i in range(len(ref_allele) - len(alt_allele)):
                                if rel_pos + len(alt_allele) < len(personalized_seq):
                                    personalized_seq.pop(rel_pos + len(alt_allele))
        
        return ''.join(personalized_seq)
    
    def __getitem__(self, key):
        """Support for genome[chrom][start:end] syntax."""
        if isinstance(key, str):
            # Return a chromosome object
            return PersonalizedChromosome(self, key)
        else:
            raise TypeError("Invalid key type")
    
    def close(self):
        """Close the reference genome."""
        self.reference_genome.close()


class PersonalizedChromosome:
    """Helper class to support genome[chrom][start:end] syntax."""
    
    def __init__(self, personalized_genome, chrom):
        self.personalized_genome = personalized_genome
        self.chrom = chrom
    
    def __getitem__(self, key):
        if isinstance(key, slice):
            return self.personalized_genome.get_sequence(self.chrom, key.start, key.stop)
        else:
            raise TypeError("Invalid key type")


def load_vcf_to_dataframe(vcf_file, sample_id):
    """
    Load VCF file into a dictionary of chromosome -> DataFrame for efficient lookup.
    
    Args:
        vcf_file: Path to VCF/BCF file
        sample_id: Sample ID to extract genotypes from
    
    Returns:
        Dictionary mapping chromosome -> DataFrame with POS, REF, ALT, and genotype columns
    """
    chromo_vcf_df = {}
    
    print(f"Loading VCF data from {vcf_file} for sample {sample_id}")
    
    with pysam.VariantFile(vcf_file, "rb") as vcf:
        if sample_id not in vcf.header.samples:
            raise ValueError(f"Sample {sample_id} not found in VCF file")
        
        for record in vcf.fetch():
            # Get genotype for the sample
            sample_gt = record.samples[sample_id]['GT']
            if sample_gt is None or sample_gt == (None, None):
                continue
            
            # Convert genotype to string format
            if "|" in str(record.samples[sample_id]):
                # Phased
                genotype_str = f"{sample_gt[0]}|{sample_gt[1]}"
            else:
                # Unphased
                genotype_str = f"{sample_gt[0]}/{sample_gt[1]}"
            
            # Create variant record
            variant_record = {
                'POS': record.pos,
                'REF': record.ref,
                'ALT': record.alts[0] if record.alts else record.ref,
                'GENOTYPE': genotype_str
            }
            
            # Add to chromosome dictionary
            chrom = record.chrom
            if chrom not in chromo_vcf_df:
                chromo_vcf_df[chrom] = []
            chromo_vcf_df[chrom].append(variant_record)
    
    # Convert lists to DataFrames for efficient filtering
    for chrom in chromo_vcf_df:
        chromo_vcf_df[chrom] = pd.DataFrame(chromo_vcf_df[chrom])
        chromo_vcf_df[chrom] = chromo_vcf_df[chrom].sort_values('POS')
    
    print(f"Loaded variants for {len(chromo_vcf_df)} chromosomes")
    return chromo_vcf_df


def get_seq_with_variants(peaks_df, genome, width, chromo_vcf_df):
    """
    Fetches sequence from genome with variant substitution, creating two haplotypes and averaging them.
    
    Args:
        peaks_df: DataFrame with chr, start, summit columns
        genome: Reference genome object
        width: Width of sequence to extract
        chromo_vcf_df: Dictionary of chromosome -> variant DataFrames
    
    Returns:
        N x L x 4 NumPy array of averaged one-hot encodings
    """
    first_vals = []
    second_vals = []
    
    for p_i, r in tqdm(peaks_df.iterrows(), total=len(peaks_df), desc="Processing sequences"):
        chromo = r['chr']
        summit = r['start'] + r['summit']
        start = summit - width // 2
        end = summit + width // 2

        # Get reference sequence
        sequence = str(genome[r['chr']][start:end])
        
        # Normalize chromosome name for VCF lookup
        if chromo.startswith('chr'):
            if chromo in ['chrX', 'chrY']:
                chromo_key = chromo
            else:
                chromo_key = int(chromo[3:])
        else:
            chromo_key = chromo

        # Get variants in this region
        c_df = chromo_vcf_df.get(chromo_key, [])
        if len(c_df) == 0:
            variants = []
        else:
            variants = c_df[(c_df['POS'] > start) & (c_df['POS'] < end)]

        # Create two haplotypes
        first_hap_seq = sequence
        second_hap_seq = sequence
        
        if len(variants) > 0:
            for i, v in variants.iterrows():
                pos = v['POS']
                ref = v['REF']
                alt = v['ALT']
                
                # Skip INDELs for now (can be extended later)
                if len(ref) > 1 or len(alt) > 1:
                    continue
                
                # Parse genotype
                genotype = v['GENOTYPE']
                
                if "|" in genotype:
                    # Phased genotype (e.g., "0|1")
                    first_hap, second_hap = genotype.split("|")
                elif "/" in genotype:
                    # Unphased genotype (e.g., "0/1")
                    first_hap, second_hap = genotype.split("/")
                else:
                    # Single allele (e.g., "1")
                    first_hap = second_hap = genotype
                
                first_hap, second_hap = int(first_hap), int(second_hap)
                
                # Convert to 0-based position relative to sequence
                pos = pos - start - 1
                
                # Validate position and reference allele
                if pos < 0 or pos + len(ref) > len(sequence):
                    continue
                    
                if sequence[pos:pos+len(ref)].upper() != ref.upper():
                    print(f"Warning: Reference mismatch at {chromo}:{pos+start+1}")
                    print(f"Expected: {ref}, Found: {sequence[pos:pos+len(ref)]}")
                    continue
                
                # Apply variants to haplotypes
                if first_hap == 1:
                    first_hap_seq = first_hap_seq[:pos] + alt + first_hap_seq[pos+len(ref):]
                if second_hap == 1:
                    second_hap_seq = second_hap_seq[:pos] + alt + second_hap_seq[pos+len(ref):]
                
                # Validate sequence lengths
                if len(first_hap_seq) != width or len(second_hap_seq) != width:
                    print(f"Warning: Sequence length mismatch after variant application")
                    continue
        
        first_vals.append(first_hap_seq.upper())
        second_vals.append(second_hap_seq.upper())

    # Convert to one-hot and average
    first_vals = one_hot.dna_to_one_hot(first_vals)
    second_vals = one_hot.dna_to_one_hot(second_vals)
    vals = np.add(first_vals, second_vals) / 2
    
    return vals


def get_seq(peaks_df, genome, width, chromo_vcf_df=None):
    """
    Fetches sequence from a given genome, with optional variant substitution.
    
    Args:
        peaks_df: DataFrame with chr, start, summit columns
        genome: Genome object (reference or personalized)
        width: Width of sequence to extract
        chromo_vcf_df: Optional VCF DataFrame for variant substitution
    """
    if chromo_vcf_df is not None:
        # Use variant-aware sequence generation
        return get_seq_with_variants(peaks_df, genome, width, chromo_vcf_df)
    else:
        # Use standard sequence generation
        vals = []
        for i, r in peaks_df.iterrows():
            sequence = str(genome[r['chr']][(r['start']+r['summit'] - width//2):(r['start'] + r['summit'] + width//2)])
            vals.append(sequence)
        return one_hot.dna_to_one_hot(vals)


def get_cts(peaks_df, bw, width):
    """
    Fetches values from a bigwig bw, given a df with minimally
    chr, start and summit columns. Summit is relative to start.
    Retrieves values of specified width centered at summit.

    "cts" = per base counts across a region
    """
    vals = []
    for i, r in peaks_df.iterrows():
        vals.append(np.nan_to_num(bw.values(r['chr'], 
                                            r['start'] + r['summit'] - width//2,
                                            r['start'] + r['summit'] + width//2)))
        
    return np.array(vals)

def get_coords(peaks_df, peaks_bool):
    """
    Fetch the co-ordinates of the regions in bed file
    returns a list of tuples with (chrom, summit)
    """
    vals = []
    for i, r in peaks_df.iterrows():
        vals.append([r['chr'], r['start']+r['summit'], "f", peaks_bool])

    return np.array(vals)

def get_seq_cts_coords(peaks_df, genome, bw, input_width, output_width, peaks_bool, chromo_vcf_df=None):
    """
    Get sequences, counts, and coordinates with optional variant substitution.
    """
    seq = get_seq(peaks_df, genome, input_width, chromo_vcf_df)
    cts = get_cts(peaks_df, bw, output_width)
    coords = get_coords(peaks_df, peaks_bool)
    return seq, cts, coords

def load_data(bed_regions, nonpeak_regions, genome_fasta, cts_bw_file, inputlen, outputlen, max_jitter, vcf_file=None, sample_id=None):
    """
    Load sequences and corresponding base resolution counts for training.
    """
    cts_bw = pyBigWig.open(cts_bw_file)
    reference_genome = pyfaidx.Fasta(genome_fasta)
    
    # Load VCF data if provided
    chromo_vcf_df = None
    if vcf_file and sample_id:
        print(f"Loading VCF data from {vcf_file} for sample {sample_id}")
        chromo_vcf_df = load_vcf_to_dataframe(vcf_file, sample_id)
    
    train_peaks_seqs = None
    train_peaks_cts = None
    train_peaks_coords = None
    train_nonpeaks_seqs = None
    train_nonpeaks_cts = None
    train_nonpeaks_coords = None

    if bed_regions is not None:
        train_peaks_seqs, train_peaks_cts, train_peaks_coords = get_seq_cts_coords(
            bed_regions, reference_genome, cts_bw, inputlen+2*max_jitter, 
            outputlen+2*max_jitter, peaks_bool=1, chromo_vcf_df=chromo_vcf_df
        )
    
    if nonpeak_regions is not None:
        train_nonpeaks_seqs, train_nonpeaks_cts, train_nonpeaks_coords = get_seq_cts_coords(
            nonpeak_regions, reference_genome, cts_bw, inputlen, outputlen, 
            peaks_bool=0, chromo_vcf_df=chromo_vcf_df
        )

    cts_bw.close()
    reference_genome.close()

    return (train_peaks_seqs, train_peaks_cts, train_peaks_coords,
            train_nonpeaks_seqs, train_nonpeaks_cts, train_nonpeaks_coords)
