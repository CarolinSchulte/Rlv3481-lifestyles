import cobra
import riptide
import csv


model = cobra.io.read_sbml_model(cur_dir+"noduleBacteria_riptide.xml")

transcriptome_tips = riptide.read_transcription_file(cur_dir+"bacteroidTips_rpkm.txt",norm=False)


with open(cur_dir+'Data/noduleBacteriaGenes.csv', newline='',encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    data = list(reader)

[nodule_genes] = data


riptide_maxfit_WT_full = riptide.maxfit_contextualize(model,transcriptome=transcriptome_tips,tasks=nodule_genes,\
                                                      frac_min=0.5,frac_max=0.95)

riptide.save_output(riptide_maxfit_WT_full,path=cur_dir+"Results/noduleBacteria")


