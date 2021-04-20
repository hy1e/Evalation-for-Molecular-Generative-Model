from utils import *
from bioactivity import jnk3_model, gsk3_model


class DistMetrics:
    def __init__(self, train_list, generated_list, on_memory=True):
        self.train_smiles_list = train_list
        self.generated_smiles_list = generated_list
        self.generated_validsmiles_list = []
        self.generated_mol_list = []
        self.jnk3_model = jnk3_model()
        self.gsk3_model = gsk3_model()

        if on_memory:
            for smiles in generated_list:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    self.generated_mol_list.append(mol)
                    self.generated_validsmiles_list.append(smiles)

    def calc_validity(self):
        return len(self.generated_validsmiles_list)/len(self.generated_smiles_list)

    def calc_novelty(self):
        num_novels = 0
        for smiles in self.generated_validsmiles_list:
            if smiles not in self.train_smiles_list:
                num_novels += 1

        return num_novels/len(self.generated_validsmiles_list)

    def calc_snn_novelty(self, r=2, th=0.4):
        fp_generated_list = [AllChem.GetMorganFingerprint(mol, r) for mol in self.generated_mol_list]
        fp_train_list = [AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smiles), r) for smiles in self.train_smiles_list]

        nov_scores = 0
        for i in range(len(fp_generated_list)):
            sim_snn = 0
            for j in range(len(fp_train_list)):
                sim = DataStructs.TanimotoSimilarity(fp_generated_list[i], fp_train_list[j])
                if sim_snn < sim:
                    sim_snn = sim

            if sim_snn < th:
                nov_scores += 1

        return nov_scores/len(fp_generated_list)

    def calc_uniqueness(self):
        num_uniques = 0
        for i in range(len(self.generated_validsmiles_list)):
            smiles = self.generated_validsmiles_list[i]
            if smiles not in self.generated_validsmiles_list[:i] + self.generated_validsmiles_list[i+1:]:
                num_uniques += 1

        return num_uniques/len(self.generated_validsmiles_list)

    def calc_diversity(self, r=2):
        fp_list = [AllChem.GetMorganFingerprint(mol, r) for mol in self.generated_mol_list]

        sim_score = 0
        for i in range(len(fp_list)):
            for j in range(i+1, len(fp_list)):
                sim_score += DataStructs.TanimotoSimilarity(fp_list[i], fp_list[j])

        return 2*sim_score/(len(fp_list)*(len(fp_list)+1))

    def calc_jnk3_act(self):
        scores = self.jnk3_model(self.generated_validsmiles_list)
        scores = np.where(scores > 0.5, 1, 0)

        return float(np.mean(scores))

    def calc_gsk3_act(self):
        scores = self.gsk3_model(self.generated_validsmiles_list)
        scores = np.where(scores > 0.5, 1, 0)

        return float(np.mean(scores))


if __name__ == '__main__':
    smiles_list = read_smilesset("Data/zinc_250k.smi")

