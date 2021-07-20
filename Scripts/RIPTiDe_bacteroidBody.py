import cobra
import riptide
import csv


model = cobra.io.read_sbml_model(cur_dir+'bacteroidBody.xml.xml')

transcriptome_body = riptide.read_transcription_file(cur_dir+'bacteroidBody_rpkm.txt',norm=False)


with open(cur_dir+'bacteroidGenes.csv', newline='',encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    data = list(reader)

[bacteroid_genes] = data


riptide_maxfit_WT_full = riptide.maxfit_contextualize(model,transcriptome=transcriptome_body,tasks=bacteroid_genes,\
                                                      frac_min=0.5,frac_max=0.95)


riptide.save_output(riptide_maxfit_WT_full,path=cur_dir+'Results/bacteroidBody')




