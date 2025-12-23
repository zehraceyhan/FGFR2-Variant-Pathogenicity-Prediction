# this predictor classifies FGFR2 protein mutations as benign or pathogenic based on sequence alignment and known variants.
# the input is an FGFR2 protein sequence in string format located in line 216.
# mutant variants can be sequenced through the "mutantvar" file and used as input for the predicor. 
from Bio import pairwise2
from blosum62_matrix import BLOSUM62_FLAT  

known_variants = {
    "V232V": "Benign",
    "L703L": "Benign",
    "T98T": "Benign",
    "L647L": "Benign",
    "A53A": "Benign",
    "V516V": "Benign",
    "T513T": "Benign",
    "T524T": "Benign",
    "Y256Y": "Benign",
    "E596E": "Benign",
    "V667V": "Benign",
    "H167H": "Benign",
    "S102S": "Benign",
    "P263P": "Benign",
    "P744P": "Benign",
    "T242T": "Benign",
    "P582P": "Benign",
    "A181A": "Benign",
    "T341T": "Benign",
    "N637N": "Benign",
    "R365R": "Benign",
    "D133D": "Benign",
    "A314A": "Benign",
    "Y237Y": "Benign",
    "A362A": "Benign",
    "T454T": "Benign",
    "S587S": "Benign",
    "H78H": "Benign",
    "G96G": "Benign",
    "L481L": "Benign",
    "H293H": "Benign",
    "R84R": "Benign",
    "F676F": "Benign",
    "F334F": "Benign",
    "L761L": "Benign",
    "A140A": "Benign",
    "S252S": "Benign",
    "Y50Y": "Benign",
    "S464S": "Benign",
    "P443P": "Benign",
    "I388I": "Benign",
    "N194N": "Benign",
    "V77V": "Benign",
    "N297N": "Benign",
    "A515A": "Benign",
    "T749T": "Benign",
    "N816N": "Benign",
    "A171A": "Benign",
    "C606C": "Benign",
    "S431S": "Benign",
    "S434S": "Benign",
    "K74K": "Benign",
    "V445V": "Benign",
    "L444L": "Benign",
    "V270V": "Benign",
    "H225H": "Benign",
    "P672P": "Benign",
    "L306L": "Benign",
    "S130S": "Benign",
    "K485K": "Benign",
    "L397L": "Benign",
    "Y308Y": "Benign",
    "T243T": "Benign",
    "E615E": "Benign",
    "C701C": "Benign",
    "C622C": "Benign",
    "Y301Y": "Benign",
    "L572L": "Benign",
    "T32T": "Benign",
    "E56E": "Benign",
    "P477P": "Benign",
    "L418L": "Benign",
    "F697F": "Benign",
    "V385V": "Benign",
    "T403T": "Benign",
    "K279K": "Benign",
    "S129S": "Benign",
    "V393V": "Benign",
    "Q752Q": "Benign",
    "V512V": "Benign",
    "T404T": "Benign",
    "V175V": "Benign",
    "R330R": "Benign",
    "T188T": "Benign",
    "H741H": "Benign",
    "G384G": "Benign",
    "S372S": "Benign",
    "S239S": "Benign",
    "D304D": "Benign",
    "A511A": "Benign",
    "E565E": "Benign",
    "I240I": "Benign",
    "A315A": "Benign",
    "Q285Q": "Benign",
    "V332V": "Benign",
    "P47P": "Benign",
    "C9C": "Benign",
    "V495V": "Benign",
    "E197E": "Benign",
    "A648T": "Pathogenic",
    "V359I": "Pathogenic",
    "A344A": "Pathogenic",
    "A344G": "Pathogenic",
    "C342R": "Pathogenic",
    "A315S": "Pathogenic",
    "W290G": "Pathogenic",
    "W290R": "Pathogenic",
    "Q289P": "Pathogenic",
    "C278Y": "Pathogenic",
    "C278F": "Pathogenic",
    "Y105C": "Pathogenic",
    "R255Q": "Pathogenic",
}


# Reference sequence = full UniProt FGFR2 sequence

REFERENCE_PROTEIN = (
    "MVSWGRFICLVVVTMATLSLARPSFSLVEDTTLEPEEPPTKYQISQPEVYVAAPGESLEVRCLLKDAAVISWTKDGVHLGPNNRTVLIGEYLQIKGATPRDSGLYACTASRTVDSETWYFMVNVTDAISSGDDEDDTDGAEDFVSENSNNKRAPYWTNTEKMEKRLHAVPAANTVKFRCPAGGNPMPTMRWLKNGKEFKQEHRIGGYKVRNQHWSLIMESVVPSDKGNYTCVVENEYGSINHTYHLDVVERSPHRPILQAGLPANASTVVGGDVEFVCKVYSDAQPHIQWIKHVEKNGSKYGPDGLPYLKVLKAAGVNTTDKEIEVLYIRNVTFEDAGEYTCLAGNSIGISFHSAWLTVLPAPGREKEITASPDYLEIAIYCIGVFLIACMVVTVILCRMKNTTKKPDFSSQPAVHKLTKRIPLRRQVTVSAESSSSMNSNTPLVRITTRLSSTADTPMLAGVSEYELPEDPKWEFPRDKLTLGKPLGEGCFGQVVMAEAVGIDKDKPKEAVTVAVKMLKDDATEKDLSDLVSEMEMMKMIGKHKNIINLLGACTQDGPLYVIVEYASKGNLREYLRARRPPGMEYSYDINRVPEEQMTFKDLVSCTYQLARGMEYLASQKCIHRDLAARNVLVTENNVMKIADFGLARDINNIDYYKKTTNGRLPVKWMAPEALFDRVYTHQSDVWSFGVLMWEIFTLGGSPYPGIPVEELFKLLKEGHRMDKPANCTNELYMMMRDCWHAVPSQRPTFKQLVEDLDRILTLTTNEEYLDLSQPLEQYSPSYPDTRSSCSSGDDSVFSPDPMPYEPCLPQYPHINGSVKT"
)





# Alignment function
def align_sequences(ref_seq, input_seq, gap_open=-10, gap_extend=-1):
    print("Starting alignment...")
    alignments = pairwise2.align.globalds(ref_seq, input_seq, BLOSUM62_FLAT, gap_open, gap_extend)
    print("Alignment finished.")
    return alignments[0]  # Best alignment


# Mutation detection
def detect_mutations(aligned_ref, aligned_input):
    mutations = []
    ref_position = 0
    

    for ref_aa, input_aa in zip(aligned_ref, aligned_input):
        if ref_aa != '-':
            ref_position += 1

        if ref_aa == '-' or input_aa == '-':
            continue  # Skip gaps


        if ref_aa != input_aa:
            mutation_id = f"{ref_aa}{ref_position}{input_aa}"
            score = BLOSUM62_FLAT.get((ref_aa, input_aa)) or BLOSUM62_FLAT.get((input_aa, ref_aa), -4)
            
            mutations.append({
                'position': ref_position,
                'ref': ref_aa,
                'alt': input_aa,
                'mutation_id': mutation_id,
                'blosum_score': score
            })
    return mutations


# Classification
def classify_mutations(mutations):
    results = []
    for m in mutations:
        mut_id = m['mutation_id']
        if mut_id in known_variants:
            label = known_variants[mut_id]
        else:
            label = "Pathogenic" if m['blosum_score'] < 0 else "Benign"

        results.append({**m, 'classification': label})
    return results


# Prediction 
def predict_from_protein(input_protein_seq):
    input_clean = input_protein_seq.replace("\n", "").replace(" ", "").strip()

    alignment = align_sequences(REFERENCE_PROTEIN, input_protein_seq)
    aligned_ref, aligned_input, score, *_ = alignment

    print("\nAligned Reference:")
    print(aligned_ref)
    print("\nAligned Input:")
    print(aligned_input)
    print(f"\nAlignment Score: {score:.2f}\n")

    mutations = detect_mutations(aligned_ref, aligned_input)

    if not mutations:
        return "No mutation detected — classified as BENIGN."

    classified = classify_mutations(mutations)

    for m in classified:
        print(
            f"Mutation: {m['mutation_id']} at position {m['position']} "
            f"| BLOSUM62 Score: {m['blosum_score']} | Classification: {m['classification']}"
        )

    if any(m['classification'] == 'Pathogenic' for m in classified):
        return "Result: At least one mutation classified as PATHOGENIC."
    else:
        return "Result: All mutations classified as BENIGN.HIYA"

# place input protein sequence below:
if __name__ == "__main__":
    input_protein = """
    # insert FGFR2 protein sequence here for prediction
    """.replace("\n", "").replace(" ", "")
    result = predict_from_protein(input_protein)
    print("\n" + result)


