#!/usr/bin/python3.6

#####!/usr/bin/env python 
# 13 Feb. 2024

# coding: utf-8

import sys
import cgi
import re, os, shutil, datetime
from collections import OrderedDict
import subprocess
import zipfile
import time
import tarfile
import glob
import random
from socket import gethostname

print ("Content-Type: text/html\n")
print ("")
#print ("Sorry, I just found this site is not working.<br>")
#print ("Ape packaged in R is not working for NJ tree reconstruction.<br>")
#print ("Could you use mirror version?<br><br>")
#print ("http://www.fish-evol.org/orthoscope/<br><br>")
#print ("Jun Inoue (26 Nov 2019)<br>")
#exit()

#time.sleep(1)

#print("hello", sys.file=stderr)
#print("hello")
#exit()
dbAddress = ""
EMAIL_SUBJECT_PREFX = '[%s]' % gethostname()
#print("EMAIL_SUBJECT_PREFX:", EMAIL_SUBJECT_PREFX)
#exit()
if "os3-373-19568.vs.sakura.ne.jp" in EMAIL_SUBJECT_PREFX:
    dbAddress  = "/var/www/html/dbb/orthoscopeDB150/"
elif "rrcs-172-254-99-39.nyc.biz.rr.com" in EMAIL_SUBJECT_PREFX:
    dbAddress   = "/dbb/orthoscopeDB150/"
elif "junINOUEp.local" in EMAIL_SUBJECT_PREFX:
    dbAddress   = "/Volumes/orthoscopeBU/dbfe/orthoscopeDB120/"
elif "yamasati" in EMAIL_SUBJECT_PREFX:
    dbAddress = "/mnt/HDD-42TB/dbb/orthoscopeDB/orthoscopeDB120/"
elif "yurai.aori.u-tokyo.ac.jp" in EMAIL_SUBJECT_PREFX:
    dbAddress = "/mnt/dbb/orthoscopeDB/orthoscopeDB150/"
elif "rx1000.site" in EMAIL_SUBJECT_PREFX:
    dbAddress = "/home1/dbb/orthoscopeDB/orthoscopeDB150/"
else:
    print("Check database address.")
    exit()

if not os.path.exists(dbAddress):
    print("Error: Database cannot be found. <br>")
    print("This server seemed to be restarted and needs database mounting. <br>")
    print("Please send an email about this error to jinoue@g.ecc.u-tokyo.ac.jp<br>")
    exit()

dirAddress = "../html/"     ##
version_orthoscope = "orthoscope"
lengthLimit_nameLine = "60"

dbLinesTMP = [
              ######## new #########
              # Sar
              ["Symbiodinium-microadriaticum_",     "000_OnlyLongestPEP_Symbiodinium_microadriaticum_gca_001939145.ASM193914v1.pep.all.fa",                "000_OnlyLongestCDS_Symbiodinium_microadriaticum_gca_001939145.ASM193914v1.cds.all.fa"],

              # Haptophyta
              ["Emiliania-huxleyi_",                "000_OnlyLongestPEP_Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.pep.all.fa",                   "000_OnlyLongestCDS_Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cds.all.fa"],
              ######## new #########

              ## Cyanobateria
              ["Prochlorococcus-marinus_",          "000_OnlyLongestPEP_Prochlorococcus_marinus_str_mit_9211.ASM1858v1.pep.all.fa",                        "000_OnlyLongestCDS_Prochlorococcus_marinus_str_mit_9211.ASM1858v1.cds.all.fa"],
              ## Actinobacteria
              ["Corynebacterium-efficiens_",        "000_OnlyLongestPEP_Corynebacterium_efficiens_ys_314.ASM1130v1.pep.all.fa",                            "000_OnlyLongestCDS_Corynebacterium_efficiens_ys_314.ASM1130v1.cds.all.fa"],
              ["Streptomyces-coelicolor_",          "000_OnlyLongestPEP_Streptomyces_coelicolor_a3_2_.ASM20383v1.pep.all.fa",                              "000_OnlyLongestCDS_Streptomyces_coelicolor_a3_2_.ASM20383v1.cds.all.fa"],
              ## Actinobacteria
              ["Methanocaldococcus-jannaschii_",    "000_OnlyLongestPEP_Methanocaldococcus_jannaschii_dsm_2661.ASM9166v1.pep.all.fa",                      "000_OnlyLongestCDS_Methanocaldococcus_jannaschii_dsm_2661.ASM9166v1.cds.all.fa"],
 
              ###Viridiplantae
              # Chlorophyta
              ["Caulerpa-lentillifera_",            "000_OnlyLongestPEP_augustus.clen170926_braker.with_hints.aa",                                         "000_OnlyLongestCDS_augustus.clen170926_braker.with_hints.fa"],

              # Charophyceae
              ["Chara-braunii_",                    "000_OnlyLongestPEP_Chara_braunii.Cbr_1.0.pep.all.fa",                                                 "000_OnlyLongestCDS_Chara_braunii.Cbr_1.0.cds.all.fa"],

              # Zygnematophyceae
              ["Ceratodon-purpureus-GG1_",          "000_CpurpureusGG1_539_v1.1.protein.fa",                                                               "000_CpurpureusGG1_539_v1.1.cds.fa"],
              ["Penium-margaritaceum_",             "000_Penium_master_pep.fa",                                                                            "000_Penium_master_CDS.fa"],

              # Marchantiophyta
              ["Marchantia-polymorpha-MpTakv61_",   "000_Marchantia-polymorpha_MpTak_v6.1r1.protein.fasta",                                                "000_Marchantia-polymorpha_MpTak_v6.1r1.cds.fasta"],
              ["Marchantia-polymorpha-MpTak1v51_",  "000_Marchantia-polymorpha_MpTak1v5.1_r2.protein.fasta",                                               "000_Marchantia-polymorpha_MpTak1v5.1_r2.cds.fasta"],
              ["Marchantia-polymorpha_",            "000_OnlyLongestPEP_Marchanta_polymorpha_v1.pep.all.fa",                                               "000_OnlyLongestCDS_Marchanta_polymorpha_v1.cds.all.fa"],


              # Anthocerotophyta
              ["Anthoceros-angustus_",              "000_Anthoceros.angustus.coding.gene.pep",                                                             "000_Anthoceros.angustus.coding.gene.cds"],
              # Bryophyta
              ["Physcomitrella-patens_",            "000_OnlyLongestPEP_Physcomitrella_patens.Phypa_V3.pep.all.fa",                                        "000_OnlyLongestCDS_Physcomitrella_patens.Phypa_V3.cds.all.fa"],


              # Lycopodiopsida 
              ["Selaginella-moellendorffii_",       "000_OnlyLongestPEP_Selaginella_moellendorffii.v1.0.pep.all.fa",                                       "000_OnlyLongestCDS_Selaginella_moellendorffii.v1.0.cds.all.fa"],
              
              # Polypodiopsida
              ["Marsilea-vestita_",                 "005_Marsilea-vestita_Mvestita_v3_protein.fa",                                                         "005_Marsilea-vestita_Mvestita_v3_cds.fa"],
              ["Azolla-filiculoides_",              "000_Azolla_filiculoides.protein.highconfidence_v1.1.fasta",                                           "000_Azolla_filiculoides.CDS.highconfidence_v1.1.fasta"],
              ["Salvinia-cucullata_",               "000_Salvinia_cucullata.protein.highconfidence_v1.2.fasta",                                            "000_Salvinia_cucullata.CDS.highconfidence_v1.2.fasta"],
              #["Alsophila-spinulosa_",              "000_OnlyLongestPEP_Als_v3.1_protein.fa",                                                              "000_OnlyLongestCDS_Als_v3.1_cds.fa"],
              ["Ceratopteris-richardii_",           "005_Ceratopteris-richardii_Crichardii_676_v2.0_protein.fa",                                           "005_Ceratopteris-richardii_Crichardii_676_v2.0_cds.fa"],

              # Acrogymnospermae
              ["Ginkgo-biloba_",                    "000_Ginkgo-biloba_proteome.gbi.csv",                                                                  "000_Ginkgo-biloba_cds.gbi.csv"],
              ["Cycas-micholitzii_",                "000_Cycas-micholitzii_proteome.cmi.csv",                                                              "000_Cycas-micholitzii_cds.cmi.csv"],
              ["Gnetum-montanum_",                  "000_Gnetum-montanum_proteome.gmo.csv",                                                                "000_Gnetum-montanum_cds.gmo.csv"],
              ["Taxus-baccata_",                    "000_Taxus-baccata_proteome.tba.csv",                                                                  "000_Taxus-baccata_cds.tba.csv"],
              ["Cryptomeria-japonica_",             "000_OnlyLongestPEP_prot_GCF_030272615.1_Sugi_1.0_rna.gbff.txt",                                       "000_OnlyLongestCDS_GCF_030272615.1_Sugi_1.0_rna.gbff.txt"],
              ["Pseudotsuga-menziesii_",            "000_Pseudotsuga-menziesii_proteome.pme.csv",                                                          "000_Pseudotsuga-menziesii_cds.pme.csv"],
              ["Picea-glauca_",                     "000_Picea-glauca_proteome.pgl.csv",                                                                   "000_Picea-glauca_cds.pgl.csv"],
              ["Pinus-sylvestris_",                 "000_Pinus-sylvestris_proteome.psy.csv",                                                               "000_Pinus-sylvestris_cds.psy.csv"],

              # Magnoliopsida
              ["Amborella-trichopoda_",             "000_OnlyLongestPEP_Amborella_trichopoda.AMTR1.0.pep.all.fa",                                          "000_OnlyLongestCDS_Amborella_trichopoda.AMTR1.0.cds.all.fa"],
              ["Nymphaea-colorata_",                "000_OnlyLongestPEP_Nymphaea_colorata.ASM883128v1.pep.all.fa",                                         "000_OnlyLongestCDS_Nymphaea_colorata.ASM883128v1.cds.all.fa"],
              ["Papaver-somniferum_",               "000_OnlyLongestPEP_Papaver_somniferum.ASM357369v1.pep.all.fa",                                        "000_OnlyLongestCDS_Papaver_somniferum.ASM357369v1.cds.all.fa"],

              # Liliopsida
              ["Dioscorea-rotundata_",              "000_OnlyLongestPEP_Dioscorea_rotundata.TDr96_F1_v2_PseudoChromosome.pep.all.fa",                      "000_OnlyLongestCDS_Dioscorea_rotundata.TDr96_F1_v2_PseudoChromosome.cds.all.fa"],
              ["Asparagus-officinalis_",            "000_OnlyLongestPEP_Asparagus_officinalis.Aspof.V1.pep.all.fa",                                        "000_OnlyLongestCDS_Asparagus_officinalis.Aspof.V1.cds.all.fa"],
              ["Musa-acuminata_",                   "000_OnlyLongestPEP_Musa_acuminata.Musa_acuminata_v2.pep.all.fa",                                      "000_OnlyLongestCDS_Musa_acuminata.Musa_acuminata_v2.cds.all.fa"],
              ["Zingiber-officinale_",              "000_OnlyLongestPEP_Zingiber-officinale__prot_GCF_018446385.1_Zo_v1.1_rna.gbff.txt",                   "000_OnlyLongestCDS_Zingiber-officinale_nucl_GCF_018446385.1_Zo_v1.1_rna.gbff.txt"],
              ["Ananas-comosus_",                   "000_OnlyLongestPEP_Ananas-comosus_nucl_prot_GCF_001540865.1_ASM154086v1_rna.gbff.txt",                "000_OnlyLongestCDS_Ananas-comosus_nucl_GCF_001540865.1_ASM154086v1_rna.gbff.txt"],

              ["Phragmites-australis_",             "000_OnlyLongestPEP_Phragmites-australis_GCF_958298935.1_lpPhrAust1.1_rna.gbff.txt",                   "000_OnlyLongestCDS_Phragmites-australis_GCF_958298935.1_lpPhrAust1.1_rna.gbff.txt"],
              ["Setaria-italica_",                  "000_OnlyLongestPEP_Setaria_italica.Setaria_italica_v2.0.pep.all.fa",                                  "000_OnlyLongestCDS_Setaria_italica.Setaria_italica_v2.0.cds.all.fa"],
              ["Sorghum-bicolor_",                  "000_OnlyLongestPEP_Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa",                                "000_OnlyLongestCDS_Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cds.all.fa"],

              ["Brachypodium-distachyon_",          "000_OnlyLongestPEP_Brachypodium_distachyon_v3.0.pep.all.fa",                                          "000_OnlyLongestCDS_Brachypodium_distachyon_v3.0.cds.all.fa"],
              ["Secale-cereale_",                   "000_OnlyLongestPEP_Secale_cereale.Rye_Lo7_2018_v1p1p1.pep.all.fa",                                    "000_OnlyLongestCDS_Secale_cereale.Rye_Lo7_2018_v1p1p1.cds.all.fa"],
              ["Hordeum-vulgare_",                  "000_OnlyLongestPEP_Hordeum_vulgare.MorexV3_pseudomolecules_assembly.pep.all.fa",                      "000_OnlyLongestCDS_Hordeum_vulgare.MorexV3_pseudomolecules_assembly.cds.all.fa"],
              ["Aegilops-tauschii_",                "000_OnlyLongestPEP_Aegilops_tauschii.Aet_v4.0.pep.all.fa",                                            "000_OnlyLongestCDS_Aegilops_tauschii.Aet_v4.0.cds.all.fa"],
              ["Triticum-spelta_",                  "000_OnlyLongestPEP_Triticum_spelta.PGSBv2.0.pep.all.fa",                                              "000_OnlyLongestCDS_Triticum_spelta.PGSBv2.0.cds.all.fa"],
              ["Triticum-aestivum_",                "000_OnlyLongestPEP_Triticum_aestivum.IWGSC.pep.all.fa",                                               "000_OnlyLongestCDS_Triticum_aestivum.IWGSC.cds.all.fa"],
              ["Leersia-perrieri_",                 "000_OnlyLongestPEP_Leersia_perrieri.Lperr_V1.4.pep.all.fa",                                           "000_OnlyLongestCDS_Leersia_perrieri.Lperr_V1.4.cds.all.fa"],
              ["Oryza-sativa-Japonica_",            "000_OnlyLongestPEP_Oryza_sativa.IRGSP-1.0.pep.all.fa",                                                "000_OnlyLongestCDS_Oryza_sativa.IRGSP-1.0.cds.all.fa"],    

              # eudicotyledons
              ["Beta-vulgaris_",                    "000_OnlyLongestPEP_Beta_vulgaris.RefBeet-1.2.2.pep.all.fa",                                           "000_OnlyLongestCDS_Beta_vulgaris.RefBeet-1.2.2.cds.all.fa"], 
              ["Actinidia-chinensis_",              "000_OnlyLongestPEP_Actinidia_chinensis.Red5_PS1_1.69.0.pep.all.fa",                                   "000_OnlyLongestCDS_Actinidia_chinensis.Red5_PS1_1.69.0.cds.all.fa"],

              # asterids
              ["Coffea-canephora_",                 "000_OnlyLongestPEP_Coffea_canephora.AUK_PRJEB4211_v1.pep.all.fa",                                     "000_OnlyLongestCDS_Coffea_canephora.AUK_PRJEB4211_v1.cds.all.fa"],    
              ["Nicotiana-attenuata_",              "000_OnlyLongestPEP_Nicotiana_attenuata.NIATTr2.pep.all.fa",                                           "000_OnlyLongestCDS_Nicotiana_attenuata.NIATTr2.cds.all.fa"],    
              ["Capsicum-annuum_",                  "000_OnlyLongestPEP_Capsicum_annuum.ASM51225v2.pep.all.fa",                                            "000_OnlyLongestCDS_Capsicum_annuum.ASM51225v2.cds.all.fa"],    
              ["Solanum-tuberosum_",                "000_OnlyLongestPEP_Solanum_tuberosum.SolTub_3.0.pep.all.fa",                                          "000_OnlyLongestCDS_Solanum_tuberosum.SolTub_3.0.cds.all.fa"],    
              ["Solanum-lycopersicum_",             "000_OnlyLongestPEP_Solanum_lycopersicum.SL3.0.pep.all.fa",                                            "000_OnlyLongestCDS_Solanum_lycopersicum.SL3.0.cds.all.fa"],    

              # rosids
              ["Vitis-vinifera_",                   "000_OnlyLongestPEP_Vitis_vinifera.12X.pep.all.fa",                                                    "000_OnlyLongestCDS_Vitis_vinifera.12X.cds.all.fa"],    
              ["Populus-trichocarpa_",              "000_OnlyLongestPEP_Populus_trichocarpa.Pop_tri_v3.pep.all.fa",                                        "000_OnlyLongestCDS_Populus_trichocarpa.Pop_tri_v3.cds.all.fa"],
              ["Cucumis-sativus_",                  "000_OnlyLongestPEP_Cucumis_sativus.ASM407v2.pep.all.fa",                                              "000_OnlyLongestCDS_Cucumis_sativus.ASM407v2.cds.all.fa"],

              ["Lupinus-angustifolius_",            "000_OnlyLongestPEP_Lupinus_angustifolius.LupAngTanjil_v1.0.pep.all.fa",                               "000_OnlyLongestCDS_Lupinus_angustifolius.LupAngTanjil_v1.0.cds.all.fa"],
              ["Trifolium-pratense_",               "000_OnlyLongestPEP_Trifolium_pratense.Trpr.pep.all.fa",                                               "000_OnlyLongestCDS_Trifolium_pratense.Trpr.cds.all.fa"],
              ["Medicago-truncatula_",              "000_OnlyLongestPEP_Medicago_truncatula.MedtrA17_4.0.pep.all.fa",                                      "000_OnlyLongestCDS_Medicago_truncatula.MedtrA17_4.0.cds.all.fa"],
              ["Phaseolus-vulgaris_",               "000_OnlyLongestPEP_Phaseolus_vulgaris.PhaVulg1_0.pep.all.fa",                                         "000_OnlyLongestCDS_Phaseolus_vulgaris.PhaVulg1_0.cds.all.fa"],
              ["Glycine-max_",                      "000_OnlyLongestPEP_Glycine_max.Glycine_max_v2.1.pep.all.fa",                                          "000_OnlyLongestCDS_Glycine_max.Glycine_max_v2.1.cds.all.fa"],
              ["Vigna-radiata_",                    "000_OnlyLongestPEP_Vigna_radiata.Vradiata_ver6.pep.all.fa",                                           "000_OnlyLongestCDS_Vigna_radiata.Vradiata_ver6.cds.all.fa"],
              ["Vigna-angularis_",                  "000_OnlyLongestPEP_Vigna_angularis.Vigan1.1.pep.all.fa",                                              "000_OnlyLongestCDS_Vigna_angularis.Vigan1.1.cds.all.fa"],

              ["Eucalyptus-grandis_",               "000_OnlyLongestPEP_Eucalyptus_grandis.Egrandis1_0.pep.all.fa",                                        "000_OnlyLongestCDS_Eucalyptus_grandis.Egrandis1_0.cds.all.fa"],
              ["Pistacia-vera_",                    "000_OnlyLongestPEP_Pistacia_vera.PisVer_v2.pep.all.fa",                                               "000_OnlyLongestCDS_Pistacia_vera.PisVer_v2.cds.all.fa"],
              ["Citrus-sinensis_",                  "000_OnlyLongestPEP_GCF_000317415.1_Csi_valencia_1.0_rna.gbff.txt",                                    "000_OnlyLongestCDS_GCF_000317415.1_Csi_valencia_1.0_rna.gbff.txt"],
              ["Theobroma-cacao_",                  "000_OnlyLongestPEP_Theobroma_cacao.Theobroma_cacao_20110822.pep.all.fa",                              "000_OnlyLongestCDS_Theobroma_cacao.Theobroma_cacao_20110822.cds.all.fa"],
              ["Gossypium-raimondii_",              "000_OnlyLongestPEP_Gossypium_raimondii.Graimondii2_0_v6.pep.all.fa",                                  "000_OnlyLongestCDS_Gossypium_raimondii.Graimondii2_0_v6.cds.all.fa"],
              ["Arabis-alpina_",                    "000_OnlyLongestPEP_Arabis_alpina.A_alpina_V4.pep.all.fa",                                             "000_OnlyLongestCDS_Arabis_alpina.A_alpina_V4.cds.all.fa"],
              ["Eutrema-salsugineum_",              "000_OnlyLongestPEP_Eutrema_salsugineum.Eutsalg1_0.pep.all.fa",                                        "000_OnlyLongestCDS_Eutrema_salsugineum.Eutsalg1_0.cds.all.fa"],
              ["Carica-papaya_",                    "000_OnlyLongestPEP_Carica-papaya_GCF_000150535.2_Papaya1.0_rna.gbff.txt",                             "000_OnlyLongestCDS_Carica-papaya_GCF_000150535.2_Papaya1.0_rna.gbff.txt"],
              ["Brassica-rapa_",                    "000_OnlyLongestPEP_Brassica_rapa.Brapa_1.0.pep.all.fa",                                               "000_OnlyLongestCDS_Brassica_rapa.Brapa_1.0.cds.all.fa"],
              ["Raphanus-sativus_",                 "000_OnlyLongestPEP_Raphanus-sativus_GCF_000801105.2_ASM80110v3_rna.gbff.txt",                         "000_OnlyLongestCDS_Raphanus-sativus_GCF_000801105.2_ASM80110v3_rna.gbff.txt"],
              ["Camelina-sativa_",                  "000_OnlyLongestPEP_Camelina_sativa.Cs.pep.all.fa",                                                    "000_OnlyLongestCDS_Camelina_sativa.Cs.cds.all.fa"], 
              ["Arabidopsis-lyrata_",               "000_OnlyLongestPEP_Arabidopsis_lyrata.v.1.0.pep.all.fa",                                              "000_OnlyLongestCDS_Arabidopsis_lyrata.v.1.0.cds.all.fa"],
              ["Arabidopsis-thaliana_",             "000_OnlyLongestPEP_Arabidopsis_thaliana.TAIR10.pep.all.fa",                                           "000_OnlyLongestCDS_Arabidopsis_thaliana.TAIR10.cds.all.fa"],


              ### Fungi
              ["Sclerotinia-sclerotiorum_",         "000_OnlyLongestPEP_Sclerotinia_sclerotiorum.ASM14694v1.pep.all.fa",                                   "000_OnlyLongestCDS_Sclerotinia_sclerotiorum.ASM14694v1.cds.all.fa"],
              ["Saccharomyces-cerevisiae_",         "000_OnlyLongestPEP_Saccharomyces_cerevisiae.R64-1-1.pep.all.fa",                                      "000_OnlyLongestCDS_Saccharomyces_cerevisiae.R64-1-1.cds.all.fa"],
              ["Saccharomyces-cerevisiae-1_",       "000_OnlyLongestPEP_Saccharomyces-cerevisiae_GCF_000146045.2_R64_rna.gbff.txt",                        "000_OnlyLongestCDS_Saccharomyces-cerevisiae_GCF_000146045.2_R64_rna.gbff.txt"],
              
              ### Ichthyosporea
              ["Sphaeroforma-arctica_",             "000_OnlyLongestPEP_S-arctica_GCF_001186125.1_Spha_arctica_JP610_V1_rna.gbff.txt",                     "000_OnlyLongestCDS_S-arctica_GCF_001186125.1_Spha_arctica_JP610_V1_rna.gbff.txt"],
              ["Capsaspora-owczarzaki_",            "000_OnlyLongestPEP_GCF_000151315.2_C_owczarzaki_V2_rna.gbff.txt",                                     "000_OnlyLongestCDS_GCF_000151315.2_C_owczarzaki_V2_rna.gbff.txt"],
              
              ## Choanoflagellida
              ["Salpingoeca-rosetta_",              "000_OnlyLongestPEP_S-rosetta_GCF_000188695.1_Proterospongia_sp_ATCC50818_rna.gbff.txt",               "000_OnlyLongestCDS_S-rosetta_GCF_000188695.1_Proterospongia_sp_ATCC50818_rna.gbff.txt"],
              ["Monosiga-brevicollis_",             "000_OnlyLongestPEP_Monosiga-brevicollis_GCF_000002865.3_V1.0_rna.gbff.txt",                           "000_OnlyLongestCDS_Monosiga-brevicollis_GCF_000002865.3_V1.0_rna.gbff.txt"],

              #### Metazoa
              ["Mnemiopsis-leidyi_",                "000_OnlyLongestPEP_Mnemiopsis_leidyi.MneLei_Aug2011.pep.all.fa",                                      "000_OnlyLongestCDS_Mnemiopsis_leidyi.MneLei_Aug2011.cds.all.fa"],
              ["Amphimedon-queenslandica_",         "000_OnlyLongestPEP_Amphimedon_queenslandica.Aqu1.pep.all.fa",                                         "000_OnlyLongestCDS_Amphimedon_queenslandica.Aqu1.cds.all.fa"],
              ["Trichoplax-adhaerens_",             "000_OnlyLongestPEP_Trichoplax_adhaerens.ASM15027v1.pep.all.fa",                                       "000_OnlyLongestCDS_Trichoplax_adhaerens.ASM15027v1.cds.all.fa"],
              
              ## Cnidaria
              ["Thelohanellus-kitauei_",            "000_OnlyLongestPEP_Thelohanellus_kitauei.ASM82789v1.pep.all.fa",                                      "000_OnlyLongestCDS_Thelohanellus_kitauei.ASM82789v1.cds.all.fa"],
              ["Hydra-vulgaris_",                   "000_OnlyLongestPEP_Hydra-vulgaris_rna.gbk.txt",                                                       "000_OnlyLongestCDS_Hydra-vulgaris_rna.gbk.txt"],
              ["Morbakka-virulenta_",               "005_Morbakka-virulenta_MOR05_r06_mRNA.longorf.pep.fa",                                                "005_Morbakka-virulenta_MOR05_r06_mRNA.longorf.cds.fa"],
              ["Aurelia-aurita-P_",                 "005_Aurelia-aurita-P_ARSv1_mRNA.longorf.pep.fa",                                                      "005_Aurelia-aurita-P_ARSv1_mRNA.longorf.cds.fa"],
              ["Aurelia-aurita-A_",                 "005_Aurelia-aurita-A_AUR21_r04_mRNA.longorf.pep.fa",                                                  "005_Aurelia-aurita-A_AUR21_r04_mRNA.longorf.cds.fa"],

              ["Dendronephthya-gigantea_",          "000_OnlyLongestPEP_Dendronephthya-gigantea_GCF_004324835.1_DenGig_1.0_rna.gbff.txt",                  "000_OnlyLongestCDS_Dendronephthya-gigantea_GCF_004324835.1_DenGig_1.0_rna.gbff.txt"],
              ["Xenia-sp_",                         "000_OnlyLongestPEP_Xenia-sp-Carnegie2017_GCF_021976095.1_XeniaSp_v1_rna.gbff.txt",                    "000_OnlyLongestCDS_Xenia-sp-Carnegie2017_GCF_021976095.1_XeniaSp_v1_rna.gbff.txt"],

              ["Nematostella-vectensis_",           "000_OnlyLongestPEP_Nematostella_vectensis.ASM20922v1.pep.all.fa",                                     "000_OnlyLongestCDS_Nematostella_vectensis.ASM20922v1.cds.all.fa"],
              ["Nematostella-vectensis-1_",         "000_OnlyLongestPEP_Nematostella-vectensis_GCF_000209225.1_ASM20922v1_rna.gbff.txt",                   "000_OnlyLongestCDS_Nematostella-vectensis_GCF_000209225.1_ASM20922v1_rna.gbff.txt"],
              ["Actinia-tenebrosa_",                "000_OnlyLongestPEP_Actinia-tenebrosa_GCF_009602425.1_ASM960242v1_rna.gbff.txt",                       "000_OnlyLongestCDS_Actinia-tenebrosa_GCF_009602425.1_ASM960242v1_rna.gbff.txt"],
              ["Exaiptasia-diaphana_",              "000_OnlyLongestPEP_Exaiptasia-diaphana_GCF_001417965.1_Aiptasia_genome_1.1_rna.gbff.txt",             "000_OnlyLongestCDS_Exaiptasia-diaphana_GCF_001417965.1_Aiptasia_genome_1.1_rna.gbff.txt"],
              ["Exaiptasia-pallida_",               "000_OnlyLongestPEP_Exaiptasia-pallida_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Exaiptasia-pallida_rna.gbk.txt"],
              ["Amplexidiscus-fenestrafer_",        "005_Amplexidiscus-fenestrafer_afen.cds.mod.longorf.pep.fa",                                           "005_Amplexidiscus-fenestrafer_afen.cds.mod.longorf.cds.fa"],
              ["Discosoma-spp_",                    "005_Discosoma-spp_dspp.cds.mod.longorf.pep.fa",                                                       "005_Discosoma-spp_dspp.cds.mod.longorf.cds.fa"],

              ["Orbicella-faveolata_",              "000_OnlyLongestPEP_Orbicella-faveolata_rna.gbk.txt",                                                  "000_OnlyLongestCDS_Orbicella-faveolata_rna.gbk.txt"],
              ["Pocillopora-verrucosa_",            "000_OnlyLongestPEP_Pocillopora-verrucosa_GCF_030620025.1_CSIRO_AGI_Pver_v1_rna.gbff.txt",             "000_OnlyLongestCDS_Pocillopora-verrucosa_GCF_030620025.1_CSIRO_AGI_Pver_v1_rna.gbff.txt"],
              ["Pocillopora-damicornis_",           "000_OnlyLongestPEP_Pocillopora-damicornis_prot_GCF_003704095.1_ASM370409v1_rna.gbff.txt",             "000_OnlyLongestCDS_Pocillopora-damicornis_nucl_GCF_003704095.1_ASM370409v1_rna.gbff.txt"],
              ["Stylophora-pistillata_",            "000_OnlyLongestPEP_Stylophora-pistillata_rna.gbk.txt",                                                "000_OnlyLongestCDS_Stylophora-pistillata_rna.gbk.txt"],
              ["Astreopora-spp_",                   "005_Astreopora-spp_astr.fa.longest-variant.pep.fa",                                                   "005_Astreopora-spp_astr.fa.longest-variant"],                    

              ["Montipora-efflorescens_",           "005_Montipora-efflorescens_meff.fa.longest-variant.pep.fa",                                           "005_Montipora-efflorescens_meff.fa.longest-variant"],  
              ["Montipora-cactus_",                 "005_Montipora-cactus_mcac.fa.longest-variant.pep.fa",                                                 "005_Montipora-cactus_mcac.fa.longest-variant"],        
              ["Acropora-tenuis_",                  "005_Acropora-tenuis_aten.fa.longest-variant.pep.fa",                                                  "005_Acropora-tenuis_aten.fa.longest-variant"],         
              ["Acropora-yongei_",                  "005_Acropora-yongei_ayon.fa.longest-variant.pep.fa",                                                  "005_Acropora-yongei_ayon.fa.longest-variant"],         
              ["Acropora-gemmifera_",               "005_Acropora-gemmifera_agem.fa.longest-variant.pep.fa",                                               "005_Acropora-gemmifera_agem.fa.longest-variant"],      
              ["Acropora-florida_",                 "005_Acropora-florida_aflo.fa.longest-variant.pep.fa",                                                 "005_Acropora-florida_aflo.fa.longest-variant"],        
              ["Acropora-awi_",                     "005_Acropora-awi_aawi.fa.longest-variant.pep.fa",                                                     "005_Acropora-awi_aawi.fa.longest-variant"],            
              ["Acropora-intermedia_",              "005_Acropora-intermedia_aint.fa.longest-variant.pep.fa",                                              "005_Acropora-intermedia_aint.fa.longest-variant"],     
              ["Acropora-acuminata_",               "005_Acropora-acuminata_aacu.fa.longest-variant.pep.fa",                                               "005_Acropora-acuminata_aacu.fa.longest-variant"],      
              ["Acropora-nasta_",                   "005_Acropora-nasta_anas.fa.longest-variant.pep.fa",                                                   "005_Acropora-nasta_anas.fa.longest-variant"],          
              ["Acropora-microphthalma_",           "005_Acropora-microphthalma_amic.fa.longest-variant.pep.fa",                                           "005_Acropora-microphthalma_amic.fa.longest-variant"],  
              ["Acropora-digitifera-1_",            "005_Acropora-digitifera_adig.fa.longest-variant.pep.fa",                                              "005_Acropora-digitifera_adig.fa.longest-variant"],   
              ["Acropora-digitifera_",              "000_OnlyLongestPEP_Acropora-digitifera_rna.gbk.txt",                                                  "000_OnlyLongestCDS_Acropora-digitifera_rna.gbk.txt"],
              ["Acropora-echinata_",                "005_Acropora-echinata_aech.fa.longest-variant.pep.fa",                                                "005_Acropora-echinata_aech.fa.longest-variant"],       
              ["Acropora-muricata_",                "005_Acropora-muricata_amur.fa.longest-variant.pep.fa",                                                "005_Acropora-muricata_amur.fa.longest-variant"],       
              ["Acropora-selago_",                  "005_Acropora-selago_asel.fa.longest-variant.pep.fa",                                                  "005_Acropora-selago_asel.fa.longest-variant"],         
              ["Acropora-hyacinthus_",              "005_Acropora-hyacinthus_ahya.fa.longest-variant.pep.fa",                                              "005_Acropora-hyacinthus_ahya.fa.longest-variant"],     
              ["Acropora-cytherea_",                "005_Acropora-cytherea_acyt.fa.longest-variant.pep.fa",                                                "005_Acropora-cytherea_acyt.fa.longest-variant"],       

              ### Spiralia
              ["Adineta-vaga_",                     "000_OnlyLongestPEP_Adineta_vaga.AMS_PRJEB1171_v1.pep.all.fa",                                         "000_OnlyLongestCDS_Adineta_vaga.AMS_PRJEB1171_v1.cds.all.fa"],
              ["Echinococcus-granulosus_",          "000_OnlyLongestPEP_Echinococcus-granulosus_GCF_000524195.1_ASM52419v1_rna.gbff.txt",                  "000_OnlyLongestCDS_Echinococcus-granulosus_GCF_000524195.1_ASM52419v1_rna.gbff.txt"],
              ["Opisthorchis-viverrini_",           "000_OnlyLongestPEP_Opisthorchis-viverrini_GCF_000715545.1_OpiViv1.0_rna.gbff.txt",                    "000_OnlyLongestCDS_Opisthorchis-viverrini_GCF_000715545.1_OpiViv1.0_rna.gbff.txt"],
              ["Schistosoma-haematobium_",          "000_OnlyLongestPEP_Schistosoma-haematobium_GCF_000699445.1_SchHae_1.0_rna.gbff.txt",                  "000_OnlyLongestCDS_Schistosoma-haematobium_GCF_000699445.1_SchHae_1.0_rna.gbff.txt"],
              ["Schistosoma-mansoni_",              "000_OnlyLongestPEP_Schistosoma_mansoni.ASM23792v2.pep.all.fa",                                        "000_OnlyLongestCDS_Schistosoma_mansoni.ASM23792v2.cds.all.fa"],
              ["Capitella-teleta_",                 "000_OnlyLongestPEP_Capitella_teleta.Capitella_teleta_v1.0.pep.all.fa",                                "000_OnlyLongestCDS_Capitella_teleta.Capitella_teleta_v1.0.cds.all.fa"],
              ["Helobdella-robusta_",               "000_OnlyLongestPEP_Helobdella_robusta.Helro1.pep.all.fa",                                             "000_OnlyLongestCDS_Helobdella_robusta.Helro1.cds.all.fa"],
              ["Notospermus-geniculatus_",          "000_OnlyLongestPEP_Notospermus-geniculatus_nge_genome_v2.0_prot.fa",                                  "000_OnlyLongestCDS_Notospermus-geniculatus_nge_genome_v2.0_nucl.fa"],
              ["Phoronis-australis_",               "000_OnlyLongestPEP_Phoronis-australis_pau_genome_v2.0_prot.fa",                                       "000_OnlyLongestCDS_Phoronis-australis_pau_genome_v2.0_nucl.fa"],
              ["Lingula-anatina_",                  "000_OnlyLongestPEP_Lingula_anatina.LinAna1.0.pep.all.fa",                                             "000_OnlyLongestCDS_Lingula_anatina.LinAna1.0.cds.all.fa"],
              ["Octopus-bimaculoides_",             "000_OnlyLongestPEP_Octopus_bimaculoides.PRJNA270931.pep.all.fa",                                      "000_OnlyLongestCDS_Octopus_bimaculoides.PRJNA270931.cds.all.fa"],
              ["Octopus-bimaculoides-1_",           "000_OnlyLongestPEP_prot_GCF_001194135.1_Octopus_bimaculoides_v2_0_rna.gbff.txt",                      "000_OnlyLongestCDS_nucl_GCF_001194135.1_Octopus_bimaculoides_v2_0_rna.gbff.txt"],
              ["Octopus-vulgaris_",                 "000_OnlyLongestPEP_out_prot_Octopus-vulgaris_GCF_006345805.1_ASM634580v1_rna.gbff.txt",               "000_OnlyLongestCDS_out_nucl_Octopus-vulgaris_GCF_006345805.1_ASM634580v1_rna.gbff.txt"],
              ["Lottia-gigantea_",                  "000_OnlyLongestPEP_Lottia_gigantea.Lotgi1.pep.all.fa",                                                "000_OnlyLongestCDS_Lottia_gigantea.Lotgi1.cds.all.fa"],
              ["Biomphalaria-glabrata_",            "000_OnlyLongestPEP_Biomphalaria-glabrata_rna.gbk.txt",                                                "000_OnlyLongestCDS_Biomphalaria-glabrata_rna.gbk.txt"],
              ["Aplysia-californica_",              "000_OnlyLongestPEP_Aplysia-californica_rna.gbk.txt",                                                  "000_OnlyLongestCDS_Aplysia-californica_rna.gbk.txt"],
              ["Pomacea-canaliculata_",             "000_OnlyLongestPEP_out_prot_Pomacea-canaliculata_GCF_003073045.1_ASM307304v1_rna.gbff.txt",           "000_OnlyLongestCDS_out_nucl_Pomacea-canaliculata_GCF_003073045.1_ASM307304v1_rna.gbff.txt"],

              #Bivalvia
              ["Mercenaria-mercenaria_",            "000_OnlyLongestPEP_Mercenaria-mercenaria_GCF_021730395.1_MADL_Memer_1_rna.gbff.txt",                  "000_OnlyLongestCDS_Mercenaria-mercenaria_GCF_021730395.1_MADL_Memer_1_rna.gbff.txt"],
              ["Mizuhopecten-yessoensis_",          "000_OnlyLongestPEP_Mizuhopecten-yessoensis_rna.gbk.txt",                                              "000_OnlyLongestCDS_Mizuhopecten-yessoensis_rna.gbk.txt"],
              ["Pecten-maximus_",                   "000_OnlyLongestPEP_Pecten-maximus_GCF_902652985.1_xPecMax1.1_rna.gbff.txt",                           "000_OnlyLongestCDS_Pecten-maximus_GCF_902652985.1_xPecMax1.1_rna.gbff.txt"],
              ["Ylistrum-balloti_",                 "000_OnlyLongestPEP_Ylistrum-balloti_GCF_031769215.1_AGI_CSIRO_Ybal_v1_rna.gbff.txt",                  "000_OnlyLongestCDS_Ylistrum-balloti_GCF_031769215.1_AGI_CSIRO_Ybal_v1_rna.gbff.txt"],
              ["Mytilus-californianus_",            "000_OnlyLongestPEP_Mytilus-californianus_GCF_021869535.1_xbMytCali1.0.p_rna.gbff.txt",                "000_OnlyLongestCDS_Mytilus-californianus_GCF_021869535.1_xbMytCali1.0.p_rna.gbff.txt"],
              ["Saccostrea-echinata_",              "000_OnlyLongestPEP_Saccostrea-echinata_GCF_033153115.1_CSIRO_AGI_Sech_v1_rna.gbff.txt",               "000_OnlyLongestCDS_Saccostrea-echinata_GCF_033153115.1_CSIRO_AGI_Sech_v1_rna.gbff.txt"],
              ["Ostrea-edulis_",                    "000_OnlyLongestPEP_Ostrea-edulis_GCF_947568905.1_xbOstEdul1.1_rna.gbff.txt",                          "000_OnlyLongestCDS_Ostrea-edulis_GCF_947568905.1_xbOstEdul1.1_rna.gbff.txt"],
              ["Crassostrea-virginica_",            "000_OnlyLongestPEP_Crassostrea-virginica_rna.gbk.txt",                                                "000_OnlyLongestCDS_Crassostrea-virginica_rna.gbk.txt"],
              ["Crassostrea-angulata_",             "000_OnlyLongestPEP_Crassostrea-angulata_GCF_025612915.1_ASM2561291v2_rna.gbff.txt",                   "000_OnlyLongestCDS_Crassostrea-angulata_GCF_025612915.1_ASM2561291v2_rna.gbff.txt"],
              ["Crassostrea-gigas_",                "000_OnlyLongestPEP_Crassostrea_gigas.oyster_v9.pep.all.fa",                                           "000_OnlyLongestCDS_Crassostrea_gigas.oyster_v9.cds.all.fa"],
              ["Pinctada-fucata_",                  "000_OnlyLongestPEP_Pinctada-fucata_pfu_aug2.0.AA1.fasta",                                             "000_OnlyLongestCDS_Pinctada-fucata_pfu_aug2.0.NT1.fasta"],
              ["Pinctada-fucata-1_",                "005_Pinctada-fucata-1_pfu_aug1.0_Pall.fasta",                                                         "005_Pinctada-fucata-1_pfu_aug1.0_Nall.fasta"],
              ["Pinctada-fucata-v41B_",             "000_OnlyLongestPEP_pfuV4.1HapB_alt_genemodels_aa.fasta",                                              "000_OnlyLongestCDS_pfuV4.1HapB_alt_genemodels_nt.fasta"],
              ["Pinctada-fucata-v41A_",             "000_OnlyLongestPEP_pfuV4.1HapA_ref_genemodels_aa.fasta",                                              "000_OnlyLongestCDS_pfuV4.1HapA_ref_genemodels_nt.fasta"],

              ## Ecdysozoa
              ["Priapulus-caudatus_",               "000_OnlyLongestPEP_Priapulus-caudatus_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Priapulus-caudatus_rna.gbk.txt"],
              ["Halicryptus-spinulosus_",           "005_PEP_Halicryptus-spinulosus_longest_orfs.pep",                                                     "005_CDS_Halicryptus-spinulosus_longest_orfs.cds"],
              ["Trichinella-spiralis_",             "000_OnlyLongestPEP_Trichinella_spiralis.Tspiralis1.pep.all.fa",                                       "000_OnlyLongestCDS_Trichinella_spiralis.Tspiralis1.cds.all.fa"],
              ["Strongyloides-ratti_",              "000_OnlyLongestPEP_Strongyloides_ratti.S_ratti_ED321_v5_0_4.pep.all.fa",                              "000_OnlyLongestCDS_Strongyloides_ratti.S_ratti_ED321_v5_0_4.cds.all.fa"],
              ["Onchocerca-volvulus_",              "000_OnlyLongestPEP_Onchocerca_volvulus.ASM49940v2.pep.all.fa",                                        "000_OnlyLongestCDS_Onchocerca_volvulus.ASM49940v2.cds.all.fa"],
              ["Loa-loa_",                          "000_OnlyLongestPEP_Loa_loa.Loa_loa_V3.pep.all.fa",                                                    "000_OnlyLongestCDS_Loa_loa.Loa_loa_V3.cds.all.fa"],
              ["Brugia-malayi_",                    "000_OnlyLongestPEP_Brugia_malayi.Bmal-4.0.pep.all.fa",                                                "000_OnlyLongestCDS_Brugia_malayi.Bmal-4.0.cds.all.fa"],
              ["Pristionchus-pacificus_",           "000_OnlyLongestPEP_Pristionchus_pacificus.P_pacificus-5.0.pep.all.fa",                                "000_OnlyLongestCDS_Pristionchus_pacificus.P_pacificus-5.0.cds.all.fa"],
              ["Caenorhabditis-japonica_",          "000_OnlyLongestPEP_Caenorhabditis_japonica.C_japonica-7.0.1.pep.all.fa",                              "000_OnlyLongestCDS_Caenorhabditis_japonica.C_japonica-7.0.1.cds.all.fa"],
              ["Caenorhabditis-brenneri_",          "000_OnlyLongestPEP_Caenorhabditis_brenneri.C_brenneri-6.0.1b.pep.all.fa",                             "000_OnlyLongestCDS_Caenorhabditis_brenneri.C_brenneri-6.0.1b.cds.all.fa"],
              ["Caenorhabditis-remanei_",           "000_OnlyLongestPEP_Caenorhabditis_remanei.C_remanei-15.0.1.pep.all.fa",                               "000_OnlyLongestCDS_Caenorhabditis_remanei.C_remanei-15.0.1.cds.all.fa"],
              ["Caenorhabditis-briggsae_",          "000_OnlyLongestPEP_Caenorhabditis_briggsae.CB4.pep.all.fa",                                           "000_OnlyLongestCDS_Caenorhabditis_briggsae.CB4.cds.all.fa"],
              ["Caenorhabditis-elegans-1_",         "000_OnlyLongestPEP_Caenorhabditis-elegans_GCF_000002985.6_WBcel235_rna.gbff.txt",                     "000_OnlyLongestCDS_Caenorhabditis-elegans_GCF_000002985.6_WBcel235_rna.gbff.txt"],
              ["Caenorhabditis-elegans_",           "000_OnlyLongestPEP_Caenorhabditis_elegans.WBcel235.pep.all.fa",                                       "000_OnlyLongestCDS_Caenorhabditis_elegans.WBcel235.cds.all.fa"],
              ["Caenorhabditis-elegans-EM58_",       "000_OnlyLongestPEP_Caenorhabditis_elegans.WBcel235.pep.all.fa",                                       "000_OnlyLongestCDS_Caenorhabditis_elegans.WBcel235.cds.all.fa"],

              ["Limulus-polyphemus_",               "000_OnlyLongestPEP_Limulus-polyphemus_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Limulus-polyphemus_rna.gbk.txt"],
              ["Sarcoptes-scabiei_",                "000_OnlyLongestPEP_Sarcoptes_scabiei.SscaA1.pep.all.fa",                                              "000_OnlyLongestCDS_Sarcoptes_scabiei.SscaA1.cds.all.fa"],
              ["Tetranychus-urticae_",              "000_OnlyLongestPEP_Tetranychus_urticae.ASM23943v1.pep.all.fa",                                        "000_OnlyLongestCDS_Tetranychus_urticae.ASM23943v1.cds.all.fa"],
              ["Ixodes-scapularis_",                "000_OnlyLongestPEP_Ixodes_scapularis.IscaW1.pep.all.fa",                                              "000_OnlyLongestCDS_Ixodes_scapularis.IscaW1.cds.all.fa"],
              ["Centruroides-sculpturatus_",        "000_OnlyLongestPEP_Centruroides-sculpturatus_rna.gbk.txt",                                            "000_OnlyLongestCDS_Centruroides-sculpturatus_rna.gbk.txt"],
              ["Parasteatoda-tepidariorum_",        "000_OnlyLongestPEP_Parasteatoda-tepidariorum_rna.gbk.txt",                                            "000_OnlyLongestCDS_Parasteatoda-tepidariorum_rna.gbk.txt"],
              ["Stegodyphus-mimosarum_",            "000_OnlyLongestPEP_Stegodyphus_mimosarum_v1.pep.all.fa",                                              "000_OnlyLongestCDS_Stegodyphus_mimosarum_v1.cds.all.fa"],
              ["Strigamia-maritima_",               "000_OnlyLongestPEP_Strigamia_maritima.Smar1.pep.all.fa",                                              "000_OnlyLongestCDS_Strigamia_maritima.Smar1.cds.all.fa"],

              # Crustacea
              ["Daphnia-magna_",                    "000_OnlyLongestPEP_Daphnia_magna.daphmag2.4.pep.all.fa",                                              "000_OnlyLongestCDS_Daphnia_magna.daphmag2.4.cds.all.fa"],
              ["Daphnia-pulex_",                    "000_OnlyLongestPEP_Daphnia_pulex.V1.0.pep.all.fa",                                                    "000_OnlyLongestCDS_Daphnia_pulex.V1.0.cds.all.fa"],
              ["Eurytemora-affinis_",               "000_OnlyLongestPEP_Eurytemora-affinis_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Eurytemora-affinis_rna.gbk.txt"],
              ["Lepeophtheirus-salmonis_",          "000_OnlyLongestPEP_Lepeophtheirus_salmonis.LSalAtl2s.pep.all.fa",                                     "000_OnlyLongestCDS_Lepeophtheirus_salmonis.LSalAtl2s.cds.all.fa"],
              ["Hyalella-azteca_",                  "000_OnlyLongestPEP_Hyalella-azteca_rna.gbk.txt",                                                      "000_OnlyLongestCDS_Hyalella-azteca_rna.gbk.txt"],
              ["Penaeus-vannamei_",                 "000_OnlyLongestPEP_Penaeus-vannamei_GCF_003789085.1_ASM378908v1_rna.gbff.txt",                        "000_OnlyLongestCDS_Penaeus-vannamei_GCF_003789085.1_ASM378908v1_rna.gbff.txt"],
              ["Penaeus-monodon_",                  "000_OnlyLongestPEP_Penaeus-monodon_GCF_015228065.1_NSTDA_Pmon_1_rna.gbff.txt",                        "000_OnlyLongestCDS_Penaeus-monodon_GCF_015228065.1_NSTDA_Pmon_1_rna.gbff.txt"],
              ["Penaeus-japonicus_",                "000_OnlyLongestPEP_Penaeus-japonicus_prot_GCF_017312705.1_Mj_TUMSAT_v1.0_rna.gbff.txt",               "000_OnlyLongestCDS_Penaeus-japonicus_nucl_GCF_017312705.1_Mj_TUMSAT_v1.0_rna.gbff.txt"],
              ["Homarus-americanus_",               "000_OnlyLongestPEP_Homarus-americanus_prot_GCF_018991925.1_GMGI_Hamer_2.0_rna.gbff.txt",              "000_OnlyLongestCDS_Homarus-americanus_nucl_GCF_018991925.1_GMGI_Hamer_2.0_rna.gbff.txt"],
              ["Portunus-trituberculatus_",         "000_OnlyLongestPEP_Portunus-trituberculatus_prot_GCF_017591435.1_ASM1759143v1_rna.gbff.txt",          "000_OnlyLongestCDS_Portunus-trituberculatus_nucl_GCF_017591435.1_ASM1759143v1_rna.gbff.txt"],

              ##### Hexapoda
              ["Pediculus-humanus_",                "000_OnlyLongestPEP_Pediculus_humanus.PhumU2.pep.all.fa",                                              "000_OnlyLongestCDS_Pediculus_humanus.PhumU2.cds.all.fa"],
              ["Zootermopsis-nevadensis_",          "000_OnlyLongestPEP_Zootermopsis_nevadensis.ZooNev1.0.pep.all.fa",                                     "000_OnlyLongestCDS_Zootermopsis_nevadensis.ZooNev1.0.cds.all.fa"],
              ["Cryptotermes-secundus_",            "000_OnlyLongestPEP_Cryptotermes-secundus_rna.gbk.txt",                                                "000_OnlyLongestCDS_Cryptotermes-secundus_rna.gbk.txt"],
              ["Rhodnius-prolixus_",                "000_OnlyLongestPEP_Rhodnius_prolixus.RproC3.pep.all.fa",                                              "000_OnlyLongestCDS_Rhodnius_prolixus.RproC3.cds.all.fa"],
              ["Acyrthosiphon-pisum_",              "000_OnlyLongestPEP_Acyrthosiphon_pisum.Acyr_2.0.pep.all.fa",                                          "000_OnlyLongestCDS_Acyrthosiphon_pisum.Acyr_2.0.cds.all.fa"],
              ["Nasonia-vitripennis_",              "000_OnlyLongestPEP_Nasonia_vitripennis.Nvit_2.1.pep.all.fa",                                          "000_OnlyLongestCDS_Nasonia_vitripennis.Nvit_2.1.cds.all.fa"],
              ["Atta-cephalotes_",                  "000_OnlyLongestPEP_Atta-cephalotes_rna.gbk.txt",                                                      "000_OnlyLongestCDS_Atta-cephalotes_rna.gbk.txt"],
              ["Solenopsis-invicta_",               "000_OnlyLongestPEP_Solenopsis-invicta_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Solenopsis-invicta_rna.gbk.txt"],
              ["Apis-mellifera_",                   "000_OnlyLongestPEP_Apis_mellifera.Amel_4.5.pep.all.fa",                                               "000_OnlyLongestCDS_Apis_mellifera.Amel_4.5.cds.all.fa"],
              ["Bombus-impatiens_",                 "000_OnlyLongestPEP_Bombus_impatiens.BIMP_2.0.pep.all.fa",                                             "000_OnlyLongestCDS_Bombus_impatiens.BIMP_2.0.cds.all.fa"],
              ["Bombus-terrestris_",                "000_OnlyLongestPEP_Bombus_terrestris.Bter_1.0.pep.all.fa",                                            "000_OnlyLongestCDS_Bombus_terrestris.Bter_1.0.cds.all.fa"],
              ["Tribolium-castaneum_",              "000_OnlyLongestPEP_Tribolium_castaneum.Tcas5.2.pep.all.fa",                                           "000_OnlyLongestCDS_Tribolium_castaneum.Tcas5.2.cds.all.fa"],
              ["Anoplophora-glabripennis_",         "000_OnlyLongestPEP_Anoplophora_glabripennis.Agla_1.0.pep.all.fa",                                     "000_OnlyLongestCDS_Anoplophora_glabripennis.Agla_1.0.cds.all.fa"],
              ["Dendroctonus-ponderosae_",          "000_OnlyLongestPEP_Dendroctonus_ponderosae.DendPond_male_1.0.pep.all.fa",                             "000_OnlyLongestCDS_Dendroctonus_ponderosae.DendPond_male_1.0.cds.all.fa"],
              ["Bombyx-mori_",                      "000_OnlyLongestPEP_Bombyx_mori.ASM15162v1.pep.all.fa",                                                "000_OnlyLongestCDS_Bombyx_mori.ASM15162v1.cds.all.fa"],
              ["Heliconius-melpomene_",             "000_OnlyLongestPEP_Heliconius_melpomene.Hmel1.pep.all.fa",                                            "000_OnlyLongestCDS_Heliconius_melpomene.Hmel1.cds.all.fa"],
              ["Danaus-plexippus_",                 "000_OnlyLongestPEP_Danaus_plexippus.Dpv3.pep.all.fa",                                                 "000_OnlyLongestCDS_Danaus_plexippus.Dpv3.cds.all.fa"],
              ["Melitaea-cinxia_",                  "000_OnlyLongestPEP_Melitaea_cinxia.MelCinx1.0.pep.all.fa",                                            "000_OnlyLongestCDS_Melitaea_cinxia.MelCinx1.0.cds.all.fa"],
              ["Belgica-antarctica_",               "000_OnlyLongestPEP_Belgica_antarctica.ASM77530v1.pep.all.fa",                                         "000_OnlyLongestCDS_Belgica_antarctica.ASM77530v1.cds.all.fa"],
              ["Anopheles-darlingi_",               "000_OnlyLongestPEP_Anopheles_darlingi.AdarC3.pep.all.fa",                                             "000_OnlyLongestCDS_Anopheles_darlingi.AdarC3.cds.all.fa"],
              ["Culex-quinquefasciatus_",           "000_OnlyLongestPEP_Culex_quinquefasciatus.CpipJ2.pep.all.fa",                                         "000_OnlyLongestCDS_Culex_quinquefasciatus.CpipJ2.cds.all.fa"],
              ["Aedes-aegypti_",                    "000_OnlyLongestPEP_Aedes_aegypti.AaegL3.pep.all.fa",                                                  "000_OnlyLongestCDS_Aedes_aegypti.AaegL3.cds.all.fa"],
              ["Mayetiola-destructor_",             "000_OnlyLongestPEP_Mayetiola_destructor.Mdes_1.0.pep.all.fa",                                         "000_OnlyLongestCDS_Mayetiola_destructor.Mdes_1.0.cds.all.fa"],
              ["Megaselia-scalaris_",               "000_OnlyLongestPEP_Megaselia_scalaris.Msca1.pep.all.fa",                                              "000_OnlyLongestCDS_Megaselia_scalaris.Msca1.cds.all.fa"],
              ["Teleopsis-dalmanni_",               "000_OnlyLongestPEP_Teleopsis_dalmanni.Tel_dalmanni_2A_v1.0.pep.all.fa",                               "000_OnlyLongestCDS_Teleopsis_dalmanni.Tel_dalmanni_2A_v1.0.cds.all.fa"],
              ["Lucilia-cuprina_",                  "000_OnlyLongestPEP_Lucilia-cuprina_rna.gbk.txt",                                                      "000_OnlyLongestCDS_Lucilia-cuprina_rna.gbk.txt"],
              ["Drosophila-grimshawi_",             "000_OnlyLongestPEP_Drosophila_grimshawi.dgri_caf1.pep.all.fa",                                        "000_OnlyLongestCDS_Drosophila_grimshawi.dgri_caf1.cds.all.fa"],
              ["Drosophila-mojavensis_",            "000_OnlyLongestPEP_Drosophila_mojavensis.dmoj_caf1.pep.all.fa",                                       "000_OnlyLongestCDS_Drosophila_mojavensis.dmoj_caf1.cds.all.fa"],
              ["Drosophila-virilis_",               "000_OnlyLongestPEP_Drosophila_virilis.dvir_caf1.pep.all.fa",                                          "000_OnlyLongestCDS_Drosophila_virilis.dvir_caf1.cds.all.fa"],
              ["Drosophila-willistoni_",            "000_OnlyLongestPEP_Drosophila_willistoni.dwil_caf1.pep.all.fa",                                       "000_OnlyLongestCDS_Drosophila_willistoni.dwil_caf1.cds.all.fa"],
              ["Drosophila-pseudoobscura_",         "000_OnlyLongestPEP_Drosophila_pseudoobscura.Dpse_3.0.pep.all.fa",                                     "000_OnlyLongestCDS_Drosophila_pseudoobscura.Dpse_3.0.cds.all.fa"],
              ["Drosophila-persimilis_",            "000_OnlyLongestPEP_Drosophila_persimilis.dper_caf1.pep.all.fa",                                       "000_OnlyLongestCDS_Drosophila_persimilis.dper_caf1.cds.all.fa"],
              ["Drosophila-ananassae_",             "000_OnlyLongestPEP_Drosophila_ananassae.dana_caf1.pep.all.fa",                                        "000_OnlyLongestCDS_Drosophila_ananassae.dana_caf1.cds.all.fa"],
              ["Drosophila-yakuba_",                "000_OnlyLongestPEP_Drosophila_yakuba.dyak_caf1.pep.all.fa",                                           "000_OnlyLongestCDS_Drosophila_yakuba.dyak_caf1.cds.all.fa"],
              ["Drosophila-erecta_",                "000_OnlyLongestPEP_Drosophila_erecta.dere_caf1.pep.all.fa",                                           "000_OnlyLongestCDS_Drosophila_erecta.dere_caf1.cds.all.fa"],
              ["Drosophila-simulans_",              "000_OnlyLongestPEP_Drosophila_simulans.ASM75419v3.pep.all.fa",                                        "000_OnlyLongestCDS_Drosophila_simulans.ASM75419v3.cds.all.fa"],
              ["Drosophila-sechellia_",             "000_OnlyLongestPEP_Drosophila_sechellia.dsec_caf1.pep.all.fa",                                        "000_OnlyLongestCDS_Drosophila_sechellia.dsec_caf1.cds.all.fa"],
              ["Drosophila-melanogaster-1_",        "000_OnlyLongestPEP_D-melanogaster_GCF_000001215.4_Release_6_plus_ISO1_MT_rna.gbff.txt",               "000_OnlyLongestCDS_D-melanogaster_GCF_000001215.4_Release_6_plus_ISO1_MT_rna.gbff.txt"],
              ["Drosophila-melanogaster_",          "000_OnlyLongestPepsTranslated_Drosophila_melanogaster.BDGP6.cds.pep.all.fa",                          "000_OnlyLongestCDStrans_Drosophila_melanogaster.BDGP6.cds.all.fa"],
              ["Drosophila-melanogaster-EM58_",      "000_OnlyLongestPEP_Drosophila_melanogaster.BDGP6.46.pep.all.fa",                                      "000_OnlyLongestCDS_Drosophila_melanogaster.BDGP6.46.cds.all.fa"],

              ### Xenacoelomorpha           
              ["Xenoturbella-bocki_",               "005_PEP_Xenoturbella-bocki_longest_orfs.pep",                                                         "005_CDS_Xenoturbella-bocki_longest_orfs.cds"],
              ["Symsagittifera-roscoffensis_",      "005_PEP_Symsagittifera-roscoffensis_longest_orfs.pep",                                                "005_CDS_Symsagittifera-roscoffensis_longest_orfs.cds"],
              #["Hofstenia-miamia_",                 "005_PEP_Hofstenia-miamia_longest_orfs.pep",                                                           "005_CDS_Hofstenia-miamia_longest_orfs.cds"],
              ["Hofstenia-miamia_",                 "000_OnlyLongestPEP_Hofstenia_miamia.HmiaM1.pep.all.fa",                                               "000_OnlyLongestCDS_Hofstenia_miamia.HmiaM1.cds.all.fa"],
              ["Isodiametra-pulchra_",              "005_PEP_Isodiametra-pulchra_longest_orfs.pep",                                                        "005_CDS_Isodiametra-pulchra_longest_orfs.cds"],
              ["Praesagittifera-naikaiensis_",      "000_PEP_Praesagittifera-naikaiensis_augustus.pna_180217c.filtered.aa.txt",                            "000_CDS_Praesagittifera-naikaiensis_augustus.pna_180217c.filtered.codingseq.txt"],

              ### Deuterostomia
              ["Saccoglossus-kowalevskii_",         "000_Saccoglossus-kowalevskii_sko.prot.txt",                                                           "000_Saccoglossus-kowalevskii_sko.cds.txt"],
              ["Saccoglossus-kowalevskii-1_",       "000_OnlyLongestPEP_Saccoglossus-kowalevskii_GCF_000003605.2_Skow_1.1_rna.gbff.txt",                   "000_OnlyLongestCDS_Saccoglossus-kowalevskii_GCF_000003605.2_Skow_1.1_rna.gbff.txt"],
              ["Ptychodera-flava_",                 "000_Ptychodera-flava_pfl.prot.txt",                                                                   "000_Ptychodera-flava_pfl.cds.txt"],
              ["Anneissia-japonica_",               "000_OnlyLongestPEP_GCF_011630105.1_ASM1163010v1_rna.gbff.txt",                                        "000_OnlyLongestCDS_GCF_011630105.1_ASM1163010v1_rna.gbff.txt"],
              ["Asterias-rubens_",                  "005_prot_GCF_902459465.1_eAstRub1.3_rna.gbff.txt",                                                    "005_nucl_GCF_902459465.1_eAstRub1.3_rna.gbff.txt"],
              ["Patiria-miniata_",                  "005_prot_GCF_015706575.1_ASM1570657v1_rna.gbff.txt",                                                  "005_nucl_GCF_015706575.1_ASM1570657v1_rna.gbff.txt"],
              ["Acanthaster-planci_",               "000_OnlyLongestPEP_Acanthaster-planci_oki-cotsv1.0.EVM2a.prot",                                       "000_OnlyLongestCDS_Acanthaster-planci_oki-cotsv1.0.EVM2a.mrna"],
              ["Acanthaster-planci-1_",             "000_OnlyLongestPEP_Acanthaster-planci_GCF_001949145.1_OKI-Apl_1.0_rna.gbff.txt",                      "000_OnlyLongestCDS_Acanthaster-planci_GCF_001949145.1_OKI-Apl_1.0_rna.gbff.txt"],
              ["Acanthaster-planci-GBR_",           "000_OnlyLongestPEP_Acanthaster-planci_gbr-cotsv1.0.EVM2a.prot",                                       "000_OnlyLongestCDS_Acanthaster-planci_gbr-cotsv1.0.EVM2a.mrna"],

              ["Lytechinus-pictus_",                "000_OnlyLongestPEP_Lytechinus-pictus_GCF_015342785.2_UCSD_Lpic_2.1_rna.gbff.txt",                     "000_OnlyLongestCDS_Lytechinus-pictus_GCF_015342785.2_UCSD_Lpic_2.1_rna.gbff.txt"],
              ["Lytechinus-variegatus_",            "000_OnlyLongestPEP_Lytechinus-variegatus_prot_GCF_018143015.1_Lvar_3.0_rna.gbff.txt",                 "000_OnlyLongestCDS_Lytechinus-variegatus_nucl_GCF_018143015.1_Lvar_3.0_rna.gbff.txt"],
              ["Strongylocentrotus-purpuratus_",    "000_OnlyLongestPEP_Strongylocentrotus_purpuratus.Spur_3.1.pep.all.fa",                                "000_OnlyLongestCDS_Strongylocentrotus_purpuratus.Spur_3.1.cds.all.fa"],
              ["Strongylocentrotus-purpuratus-1_",  "000_OnlyLongestPEP_S-purpuratus_GCF_000002235.4_Spur_4.2_rna.gbff.txt",                               "000_OnlyLongestCDS_S-purpuratus_GCF_000002235.4_Spur_4.2_rna.gbff.txt"],
              ["Strongylocentrotus-purpuratus-S50_","000_OnlyLongestPEP_Strongylocentrotus-purpuratus-Spur50_GCF_000002235.5_Spur_5.0_rna.gbff.txt",       "000_OnlyLongestCDS_Strongylocentrotus-purpuratus-Spur50_GCF_000002235.5_Spur_5.0_rna.gbff.txt"],

              ["Asymmetron-lucayanum-CH_",          "005_Asymmetron-lucayanum-CH_cdhit-prot98.fs",                                                         "005_Asymmetron-lucayanum-CH_cdhit-nucl98.fs"],
              ["Asymmetron-lucayanum_",             "005_Asymmetron-lucayanum_longest_orfs.pep",                                                           "005_Asymmetron-lucayanum_longest_orfs.cds"],

              #["Branchiostoma-belcheri_",           "000_OnlyLongestPEP_Branchiostoma.belcheri_v18h27.r3_Translated.fa",                                   "000_OnlyLongestCDS_Branchiostoma.belcheri_v18h27.r3_ref_cds.fa"],
              #["Branchiostoma-belcheri-1_",         "000_OnlyLongestPEP_B-belcheri_GCF_001625305.1_Haploidv18h27_rna.gbff.txt",                            "000_OnlyLongestCDS_B-belcheri_GCF_001625305.1_Haploidv18h27_rna.gbff.txt"],
              #["Branchiostoma-belcheri-1_",         "000_OnlyLongestPEP_GCF_001625305.1_Haploidv18h27_rna.gbff.txt",                                       "000_OnlyLongestCDS_GCF_001625305.1_Haploidv18h27_rna.gbff.txt"],
              ["Branchiostoma-belcheri-RS2_",       "000_OnlyLongestPEP_Branchiostoma-belcheri-rs2_GCF_001625305.1_Haploidv18h27_rna.gbff.txt",            "000_OnlyLongestCDS_Branchiostoma-belcheri-rs2_GCF_001625305.1_Haploidv18h27_rna.gbff.txt"],

              ["Branchiostoma-japonicum-1_",        "000_OnlyLongestPEP_Branchiostoma-japonicum_aa69938.txt",                                              "000_OnlyLongestCDS_Branchiostoma-japonicum_cds69938.txt"],
              ["Branchiostoma-japonicum_",          "000_OnlyLongestPEP_Branchiostoma-japonicum_selected_amphioxus.polished.aligned10917_expression.aa",   "000_OnlyLongestCDS_Branchiostoma-japonicum_selected_amphioxus.polished.aligned10917_expression.fasta"],
              #["Branchiostoma-japonicum_",          "000_OnlyLongestPEP_Branchiostoma-japonicum_aa69938.txt",                                              "000_OnlyLongestCDS_Branchiostoma-japonicum_cds69938.txt"],

              ["Branchiostoma-lanceolatum_",        "000_OnlyLongestPEP_Bla_annot_v4_best_Aac.fa",                                                         "000_OnlyLongestCDS_Bla_annot_v4_best_Cds.fa"],
              ["Branchiostoma-lanceolatum-E_",      "000_OnlyLongestPEP_Branchiostoma_lanceolatum.BraLan2.pep.all.fa",                                     "000_OnlyLongestCDS_Branchiostoma_lanceolatum.BraLan2.cds.all.fa"],
              #["Branchiostoma-floridae_",           "000_OnlyLongestPEP_Branchiostoma-floridae_proteins.Brafl1_21954.fasta",                               "000_OnlyLongestCDS_Branchiostoma-floridae_transcripts.Brafl1_21954.fasta"],
              ["Branchiostoma-floridae-1_",         "000_OnlyLongestPEP_B-floridae_GCF_000003815.1_Version_2_rna.gbff.txt",                                "000_OnlyLongestCDS_B-floridae_GCF_000003815.1_Version_2_rna.gbff.txt"],
              #["Branchiostoma-floridae-rs_",        "000_OnlyLongestPEP_GCF_000003815.2_Bfl_VNyyK_rna.gbff.txt",                                           "000_OnlyLongestCDS_GCF_000003815.2_Bfl_VNyyK_rna.gbff.txt"],
              ["Branchiostoma-floridae-RS2_",       "000_OnlyLongestPEP_Branchiostoma-floridae-rs2_GCF_000003815.2_Bfl_VNyyK_rna.gbff.txt",                "000_OnlyLongestCDS_Branchiostoma-floridae-rs2_GCF_000003815.2_Bfl_VNyyK_rna.gbff.txt"],

              ["Fritillaria-haplostoma_",           "005_Hososaizuchiboya_prot.txt",                                                                       "005_Hososaizuchiboya_nucl.txt"],   
              ["Oikopleura-rufescens_",             "005_Maruotamaboya_prot.txt",                                                                          "005_Maruotamaboya_nucl.txt"],     
              ["Oikopleura-longicauda_",            "005_Onagaotamaboya_prot.txt",                                                                         "005_Onagaotamaboya_nucl.txt"],    
              ["Oikopleura-dioica_",                "000_OnlyLongestPEP_Oikopleura_peptides_reference_v1.0a.fa.txt",                                       "000_OnlyLongestCDS_Oikopleura_transcripts_reference_v1.0.fa.txt"],

              ["Dilioletta-gegenbauri_",            "005_Ooumitaru_prot.txt",                                                                              "005_Ooumitaru_nucl.txt"],    
              ["Doliolum-denticulatum_",            "005_Doliolum-denticulatum_prot1.txt",                                                                 "005_Doliolum-denticulatum_nucl1.txt"],
              ["Doliolum-denticulatum-L_",          "005_Doliolum-denticulatum-L_prot.txt",                                                                "005_Doliolum-denticulatum-L_nucl.txt"],
              ["Doliolum-nationalis_",              "005_Doliolum-nationalis_longest_orfs.pep",                                                            "005_Doliolum-nationalis_longest_orfs.cds"],
              ["Doliolum-nationalis-I_",            "005_Himeumitaru_prot.txt",                                                                            "005_Himeumitaru_nucl.txt"],

              ["Thalia-orientalis_",                "005_Himesaruba_prot.txt",                                                                             "005_Himesaruba_nucl.txt"],        
              ["Salpa-thompsoni-t_",                "005_tSalTho2.1_Genome_masked_protein.fa",                                                             "005_tSalTho2.1_Genome_masked_cds.fa"],
              ["Salpa-thompsoni_",                  "005_Salpa-thompsoni_longest_orfs.pep",                                                                "005_Salpa-thompsoni_longest_orfs.cds"],
              ["Salpa-fusiformis_",                 "005_Salpa_fusiformis_SRR6326577_longest_orfs.pep",                                                    "005_Salpa_fusiformis_SRR6326577_longest_orfs.cds"],

              ["Phallusia-mammillata_",             "005_Phallusia-mammillata_longest_orfs.pep",                                                           "005_Phallusia-mammillata_longest_orfs.cds"],
              ["Ciona-savignyi_",                   "000_OnlyLongestPEP_Ciona_savignyi.CSAV2.0.pep.all.fa",                                                "000_OnlyLongestCDS_Ciona_savignyi.CSAV2.0.cds.all.fa"],
              ["Ciona-savignyi-1_",                 "000_OnlyLongestPEP_Ciona_savignyi.CSAV2.0.81_pep.renamed.fa",                                         "000_OnlyLongestCDS_Ciona_savignyi.CSAV2.0.cdna.ncdna_renamedGoodNomenclature.fa"],
              ["Ciona-intestinalis_",               "000_OnlyLongestPEP_Ciona_intestinalis.KH.pep.all.fa",                                                 "000_OnlyLongestCDS_Ciona_intestinalis.KH.cds.all.fa"],
              ["Ciona-intestinalis-1_",             "000_OnlyLongestPEP_Ciona-intestinalis_GCF_000224145.3_KH_rna.gbff.txt",                               "000_OnlyLongestCDS_Ciona-intestinalis_GCF_000224145.3_KH_rna.gbff.txt"],
              ["Clavelina-lepadiformis_",           "005_Clavelina-lepadiformis_longest_orfs.pep",                                                         "005_Clavelina-lepadiformislongest_orfs.cds"],
              ["Cystodytes-dellechiajei_",          "005_Cystodytes-dellechiajei_longest_orfs.pep",                                                        "005_Cystodytes-dellechiajei_longest_orfs.cds"],

              ["Bostrichobranchus-pilularis_",      "005_Bostrichobranchus-pilularis_longest_orfs.pep",                                                    "005_Bostrichobranchus-pilularis_longest_orfs.cds"],
              ["Molgula-occidentalis_",             "005_moxi-maker2.all.maker.proteins_goodNomenclature_filtered.aa",                                     "005_moxi-maker2.all.maker.transcripts_goodNomenclature_filtered.fasta"],
              ["Molgula-occidentalis-1_",           "005_Molgula_occidentalis-1_longest_orfs.pep",                                                         "005_Molgula_occidentalis-1_longest_orfs.cds"],
              ["Molgula-manhattensis_",             "005_Molgula-manhattensis_longest_orfs.pep",                                                           "005_Molgula-manhattensis_longest_orfs.cds"],
              ["Molgula-oculata_",                  "005_PEP_M_oculata_maker2.all9_proteins_goodNomenclature_filtered.aa",                                 "005_CDS_M_oculata_maker2.all9_transcripts_goodNomenclature_filtered.fasta"],

              ["Styela-plicata_",                   "005_Styela-plicata_longest_orfs.pep",                                                                 "005_Styela-plicata_longest_orfs.cds"],
              ["Styela-clava_",                     "000_OnlyLongestPEP_Styela-clava_GCF_013122585.1_ASM1312258v2_rna.gbff.txt",                           "000_OnlyLongestCDS_Styela-clava_GCF_013122585.1_ASM1312258v2_rna.gbff.txt"],
              ["Botrylloides-leachii_",             "005_PEP_Boleac_proteins_v4.fasta",                                                                    "005_CDS_Boleac_transcripts_v5.fasta"],
              ["Botryllus-schlosseri_",             "000_PEP_Botryllus-schlosseri_start_stop_transcripts30.fa.txt",                                        "000_CDS_Botryllus-schlosseri_start_stop_transcripts30.fa.txt"],
              ["Dendrodoa-grossularia_",            "005_Dendrodoa-grossularia_longest_orfs.pep",                                                          "005_Dendrodoa-grossularia_longest_orfs.cds"],
              ["Polyandrocarpa-anguinea_",          "005_Polyandrocarpa-anguinea_longest_orfs.pep",                                                        "005_Polyandrocarpa-anguinea_longest_orfs.cds"],
              ["Microcosmus-squamiger_",            "005_Microcosmus-squamiger_longest_orfs.pep",                                                          "005_Microcosmus-squamiger__longest_orfs.cds"],
              ["Halocynthia-aurantium_",            "000_OnlyLongestPEP_Halocynthia-aurantium_Haaura.MTP2014.proteins_grouped.aa",                         "000_OnlyLongestCDS_Halocynthia-aurantium_Haaura.MTP2014_transcripts_grouped.fasta"],
              ["Halocynthia-roretzi_",              "000_OnlyLongestPEP_Halocynthia-roretzi_Harore.MTP2014.protein_2018.aa",                               "000_OnlyLongestCDS_Halocynthia-roretzi_Harore.MTP2014.transcript_2018.fasta"],

              ### Vertebrates
              ["Eptatretus-burgeri_",               "000_OnlyLongestPEP_Eptatretus_burgeri.Eburgeri_3.2.pep.all.fa",                                       "000_OnlyLongestCDS_Eptatretus_burgeri.Eburgeri_3.2.cds.all.fa"],
              ["Petromyzon-marinus_",               "000_Petromyzon-marinus_PMZ_v3.1_proteins.fa",                                                         "000_Petromyzon-marinus_PMZ_v3.1_transcripts.fa"],

              ["Callorhinchus-milii_",              "000_OnlyLongestPEP_GCF_000165045.1_Callorhinchus_milii-6.1.3_rna.gbff.txt",                           "000_OnlyLongestCDS_GCF_000165045.1_Callorhinchus_milii-6.1.3_rna.gbff.txt"],
              ["Callorhinchus-milii-E102_",         "000_OnlyLongestPEP_Callorhinchus_milii.Callorhinchus_milii-6.1.3.pep.all.fa",                         "000_OnlyLongestCDS_Callorhinchus_milii.Callorhinchus_milii-6.1.3.cds.all.fa"],
              ["Leucoraja-erinacea_",               "000_OnlyLongestPEP_GCF_028641065.1_Leri_hhj_1_rna.gbff.txt",                                          "000_OnlyLongestCDS_GCF_028641065.1_Leri_hhj_1_rna.gbff.txt"],
              ["Amblyraja-radiata-N_",              "000_OnlyLongestPEP_Amblyraja-radiata_GCF_010909765.1_sAmbRad1.pri_rna.gbff.txt",                      "000_OnlyLongestCDS_Amblyraja-radiata_GCF_010909765.1_sAmbRad1.pri_rna.gbff.txt"],
              ["Pristis-pectinata_",                "000_OnlyLongestPEP_GCF_009764475.1_sPriPec2.1.pri_rna.gbff.txt",                                      "000_OnlyLongestCDS_GCF_009764475.1_sPriPec2.1.pri_rna.gbff.txt"],
              ["Mobula-hypostoma_",                 "000_OnlyLongestPEP_Mobula-hypostoma_GCF_963921235.1_sMobHyp1.1_rna.gbff.txt",                         "000_OnlyLongestCDS_Mobula-hypostoma_GCF_963921235.1_sMobHyp1.1_rna.gbff.txt"],
              ["Hypanus-sabinus_",                  "000_OnlyLongestPEP_GCF_030144855.1_sHypSab1.hap1_rna.gbff.txt",                                       "000_OnlyLongestCDS_GCF_030144855.1_sHypSab1.hap1_rna.gbff.txt"],
              ["Hemiscyllium-ocellatum_",           "000_OnlyLongestPEP_Hemiscyllium-ocellatum_GCF_020745735.1_sHemOce1.pat.X.cur._rna.gbff.txt",          "000_OnlyLongestCDS_Hemiscyllium-ocellatum_GCF_020745735.1_sHemOce1.pat.X.cur._rna.gbff.txt"],
              ["Chiloscyllium-punctatum_",          "005_Cpunctatum_v1.0.pep.faa",                                                                         "005_Cpunctatum_v1.0.cds.nuc.fna"],
              ["Chiloscyllium-plagiosum_",          "000_OnlyLongestPEP_Chiloscyllium-plagiosum_GCF_004010195.1_ASM401019v2_rna.gbff.txt",                 "000_OnlyLongestCDS_Chiloscyllium-plagiosum_GCF_004010195.1_ASM401019v2_rna.gbff.txt"],
              ["Stegostoma-fasciatum_",             "000_OnlyLongestPEP_Stegostoma-fasciatum_GCF_022316705.1_sSteFas1.1_rna.gbff.txt",                     "000_OnlyLongestCDS_Stegostoma-fasciatum_GCF_022316705.1_sSteFas1.1_rna.gbff.txt"],
              ["Stegostoma-tigrinum_",              "000_OnlyLongestPEP_GCF_030684315.1_sSteTig4.hap1_rna.gbff.txt",                                       "000_OnlyLongestCDS_GCF_030684315.1_sSteTig4.hap1_rna.gbff.txt"],
              ["Rhincodon-typus_",                  "000_OnlyLongestPEP_Rhincodon-typus_rna.gbk.txt",                                                      "000_OnlyLongestCDS_Rhincodon-typus_rna.gbk.txt"],
              ["Rhincodon-typus-R_",                "005_Rtypus_kobe_v1.0.pep.faa",                                                                        "005_Rtypus_kobe_v1.0.cds.nuc.fna"],
              ["Carcharodon-carcharias_",           "000_OnlyLongestPEP_GCF_017639515.1_sCarCar2.pri_rna.gbff.txt",                                        "000_OnlyLongestCDS_GCF_017639515.1_sCarCar2.pri_rna.gbff.txt"],
              ["Scyliorhinus-torazame_",            "005_Storazame_v1.0.pep.faa",                                                                          "005_Storazame_v1.0.cds.nuc.fna"],
              ["Scyliorhinus-canicula_",            "000_OnlyLongestPEP_GCF_902713615.1_sScyCan1.1_rna.gbff.txt",                                          "000_OnlyLongestCDS_GCF_902713615.1_sScyCan1.1_rna.gbff.txt"],

              ### Actinopterygii
              ["Erpetoichthys-calabaricus_",        "000_OnlyLongestPEP_prot_Erpetoichthys-calabaricus_GCF_900747795.1_fErpCal1.1_rna.gbff.txt",           "000_OnlyLongestCDS_nucl_Erpetoichthys-calabaricus_GCF_900747795.1_fErpCal1.1_rna.gbff.txt"],
              ["Polypterus-senegalus_",             "000_OnlyLongestPEP_Polypterus-senegalus_prot_GCF_016835505.1_ASM1683550v1_rna.gbff.txt",              "000_OnlyLongestCDS_Polypterus-senegalus_nucl_GCF_016835505.1_ASM1683550v1_rna.gbff.txt"],

              ["Polyodon-spathula_",                "000_OnlyLongestPEP_Polyodon-spathula_prot_GCF_017654505.1_ASM1765450v1_rna.gbff.txt",                 "000_OnlyLongestCDS_Polyodon-spathula_nucl_GCF_017654505.1_ASM1765450v1_rna.gbff.txt"],
              ["Acipenser-ruthenus_",               "000_OnlyLongestPEP_prot_A-ruthenus_GCF_010645085.1_ASM1064508v1_rna.gbff.txt",                        "000_OnlyLongestCDS_nucl_A-ruthenus_GCF_010645085.1_ASM1064508v1_rna.gbff.txt"],

              ["Lepisosteus-oculatus_",             "000_OnlyLongestPEP_Lepisosteus_oculatus.LepOcu1.pep.all.fa",                                          "000_OnlyLongestCDS_Lepisosteus_oculatus.LepOcu1.cds.all.fa"],
              ["Lepisosteus-oculatus-1_",           "000_OnlyLongestPEP_Lepisosteus-oculatus_GCF_000242695.1_LepOcu1_rna.gbff.txt",                        "000_OnlyLongestCDS_Lepisosteus-oculatus_GCF_000242695.1_LepOcu1_rna.gbff.txt"],

              ["Amia-calva_",                       "000_Amia-calva_AMC_proteins.fasta",                                                                   "000_Amia-calva_AMC_transcripts.fasta"],

              ["Scleropages-formosus_",             "000_OnlyLongestPEP_Scleropages-formosus_rna.gbk.txt",                                                 "000_OnlyLongestCDS_Scleropages-formosus_rna.gbk.txt"],
              ["Paramormyrops-kingsleyae_",         "000_OnlyLongestPEP_Paramormyrops-kingsleyae_rna.gbk.txt",                                             "000_OnlyLongestCDS_Paramormyrops-kingsleyae_rna.gbk.txt"],
              ["Megalops-cyprinoides_",             "000_OnlyLongestPEP_GCF_013368585.1_fMegCyp1.pri_rna.gbff.txt",                                        "000_OnlyLongestCDS_GCF_013368585.1_fMegCyp1.pri_rna.gbff.txt"],
              ["Conger-conger_",                    "000_OnlyLongestPEP_prot_GCF_963514075.1_fConCon1.1_rna.gbff.txt",                                     "000_OnlyLongestCDS_nucl_GCF_963514075.1_fConCon1.1_rna.gbff.txt"],
              ["Anguilla-anguilla_",                "000_OnlyLongestPEP_prot_GCF_013347855.1_fAngAng1.pri_rna.gbff.txt",                                   "000_OnlyLongestCDS_nucl_GCF_013347855.1_fAngAng1.pri_rna.gbff.txt"],

              ["Denticeps-clupeoides_",             "000_OnlyLongestPEP_prot_Denticeps-clupeoides_GCF_900700375.1_fDenClu1.1_rna.gbff.txt",                "000_OnlyLongestCDS_nucl_Denticeps-clupeoides_GCF_900700375.1_fDenClu1.1_rna.gbff.txt"],
              ["Clupea-harengus_",                  "000_OnlyLongestPEP_Clupea-harengus_rna.gbk.txt",                                                      "000_OnlyLongestCDS_Clupea-harengus_rna.gbk.txt"],
              ["Danio-rerio_",                      "000_OnlyLongestPEP_Danio_rerio.GRCz10.pep.all.fa",                                                    "000_OnlyLongestCDS_Danio_rerio.GRCz10.cds.all.fa"],
              ["Danio-rerio-1_",                    "000_OnlyLongestPEP_Danio-rerio_zebrafish.1.rna.gbff.txt",                                             "000_OnlyLongestCDS_Danio-rerio_zebrafish.1.rna.gbff.txt"],
              ["Carassius-auratus_",                "000_OnlyLongestPEP_prot_Carassius-auratus_GCF_003368295.1_ASM336829v1_rna.gbff.txt",                  "000_OnlyLongestCDS_nucl_Carassius-auratus_GCF_003368295.1_ASM336829v1_rna.gbff.txt"],
              ["Sinocyclocheilus-grahami_",         "000_OnlyLongestPEP_Sinocyclocheilus-grahami_rna.gbk.txt",                                             "000_OnlyLongestCDS_Sinocyclocheilus-grahami_rna.gbk.txt"],
              ["Cyprinus-carpio_",                  "000_OnlyLongestPEP_Cyprinus-carpio_rna.gbk.txt",                                                      "000_OnlyLongestCDS_Cyprinus-carpio_rna.gbk.txt"],
              ["Electrophorus-electricus_",         "000_OnlyLongestPEP_prot_Electrophorus-electricus_GCF_003665695.1_Ee_SOAP_WITH_SSPACE_rna.gbff.txt",   "000_OnlyLongestCDS_nucl_Electrophorus-electricus_GCF_003665695.1_Ee_SOAP_WITH_SSPACE_rna.gbff.txt"],
              ["Tachysurus-fulvidraco_",            "000_OnlyLongestPEP_prot_Tachysurus-fulvidraco_GCF_003724035.1_ASM372403v1_rna.gbff.txt",              "000_OnlyLongestCDS_nucl_Tachysurus-fulvidraco_GCF_003724035.1_ASM372403v1_rna.gbff.txt"],
              ["Pangasianodon-hypophthalmus_",      "000_OnlyLongestPEP_prot_Pangasianodon-hypophthalmus_GCF_003671635.1_VN_pangasius_rna.gbff.txt",       "000_OnlyLongestCDS_nucl_Pangasianodon-hypophthalmus_GCF_003671635.1_VN_pangasius_rna.gbff.txt"],
              ["Ictalurus-punctatus_",              "000_OnlyLongestPEP_Ictalurus-punctatus_rna.gbk.txt",                                                  "000_OnlyLongestCDS_Ictalurus-punctatus_rna.gbk.txt"],
              ["Astyanax-mexicanus_",               "000_OnlyLongestPEP_Astyanax_mexicanus.AstMex102.pep.all.fa",                                          "000_OnlyLongestCDS_Astyanax_mexicanus.AstMex102.cds.all.fa"],
              ["Pygocentrus-nattereri_",            "000_OnlyLongestPEP_Pygocentrus-nattereri_rna.gbk.txt",                                                "000_OnlyLongestCDS_Pygocentrus-nattereri_rna.gbk.txt"],
              ["Esox-lucius_",                      "000_OnlyLongestPEP_Esox-lucius_rna.gbk.txt",                                                          "000_OnlyLongestCDS_Esox-lucius_rna.gbk.txt"],
              ["Hucho-hucho_",                      "000_OnlyLongestPEP_Hucho_hucho.ASM331708v1.pep.all.fa",                                               "000_OnlyLongestCDS_Hucho_hucho.ASM331708v1.cds.all.fa"],

              ["Coregonus-clupeaformis_",           "000_OnlyLongestPEP_Coregonus-clupeaformis_GCF_020615455.1_ASM2061545v1_rna.gbff.txt",                 "000_OnlyLongestCDS_Coregonus-clupeaformis_GCF_020615455.1_ASM2061545v1_rna.gbff.txt"],
              ["Salmo-trutta_",                     "000_OnlyLongestPEP_Salmo-trutta_GCF_901001165.1_fSalTru1.1_rna.gbff.txt",                             "000_OnlyLongestCDS_Salmo-trutta_GCF_901001165.1_fSalTru1.1_rna.gbff.txt"],
              ["Salmo-salar_",                      "000_OnlyLongestPEP_Salmo-salar_rna.gbk.txt",                                                          "000_OnlyLongestCDS_Salmo-salar_rna.gbk.txt"],
              ["Salvelinus-namaycush_",             "000_OnlyLongestPEP_Salvelinus-namaycush_GCF_016432855.1_SaNama_1.0_rna.gbff.txt",                     "000_OnlyLongestCDS_Salvelinus-namaycush_GCF_016432855.1_SaNama_1.0_rna.gbff.txt"],
              ["Salvelinus-alpinus_",               "000_OnlyLongestPEP_Salvelinus-alpinus_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Salvelinus-alpinus_rna.gbk.txt"],
              ["Oncorhynchus-mykiss_",              "000_OnlyLongestPEP_Oncorhynchus-mykiss_rna.gbk.txt",                                                  "000_OnlyLongestCDS_Oncorhynchus-mykiss_rna.gbk.txt"],
              ["Oncorhynchus-tshawytscha_",         "000_OnlyLongestPEP_Oncorhynchus-tshawytscha_rna.gbk.txt",                                             "000_OnlyLongestCDS_Oncorhynchus-tshawytscha_rna.gbk.txt"],
              ["Oncorhynchus-kisutch_",             "000_OnlyLongestPEP_Oncorhynchus-kisutch_rna.gbk.txt",                                                 "000_OnlyLongestCDS_Oncorhynchus-kisutch_rna.gbk.txt"],
              ["Oncorhynchus-nerka_",               "000_OnlyLongestPEP_Oncorhynchus-nerka_GCF_006149115.2_Oner_1.1_rna.gbff.txt",                         "000_OnlyLongestCDS_Oncorhynchus-nerka_GCF_006149115.2_Oner_1.1_rna.gbff.txt"],
              ["Oncorhynchus-gorbuscha_",           "000_OnlyLongestPEP_Oncorhynchus-gorbuscha_GCF_021184085.1_OgorEven_v1.0_rna.gbff.txt",                "000_OnlyLongestCDS_Oncorhynchus-gorbuscha_GCF_021184085.1_OgorEven_v1.0_rna.gbff.txt"],
              ["Oncorhynchus-keta_",                "000_OnlyLongestPEP_Oncorhynchus-keta_GCF_023373465.1_Oket_V2_rna.gbff.txt",                           "000_OnlyLongestCDS_Oncorhynchus-keta_GCF_023373465.1_Oket_V2_rna.gbff.txt"],
              ["Gadus-morhua_",                     "000_OnlyLongestPEP_Gadus_morhua.gadMor1.pep.all.fa",                                                  "000_OnlyLongestCDS_Gadus_morhua.gadMor1.cds.all.fa"],
              ["Labrus-bergylta_",                  "000_OnlyLongestPEP_Labrus-bergylta_rna.gbk.txt",                                                      "000_OnlyLongestCDS_Labrus-bergylta_rna.gbk.txt"],
              ["Sander-lucioperca_",                "000_OnlyLongestPEP_Sander_lucioperca.SLUC_FBN_1.pep.all.fa",                                          "000_OnlyLongestCDS_Sander_lucioperca.SLUC_FBN_1.cds.all.fa"],
              ["Cyclopterus-lumpus_",               "000_OnlyLongestPEP_Cyclopterus_lumpus.fCycLum1.pri.pep.all.fa",                                       "000_OnlyLongestCDS_Cyclopterus_lumpus.fCycLum1.pri.cds.all.fa"],
              ["Gasterosteus-aculeatus_",           "000_OnlyLongestPEP_Gasterosteus_aculeatus.BROADS1.pep.all.fa",                                        "000_OnlyLongestCDS_Gasterosteus_aculeatus.BROADS1.cds.all.fa"],
              ["Cottoperca-gobio_",                 "000_OnlyLongestPEP_Cottoperca_gobio.fCotGob3.1.pep.all.fa",                                           "000_OnlyLongestCDS_Cottoperca_gobio.fCotGob3.1.cds.all.fa"],
              ["Notothenia-coriiceps_",             "000_OnlyLongestPEP_Notothenia-coriiceps_rna.gbk.txt",                                                 "000_OnlyLongestCDS_Notothenia-coriiceps_rna.gbk.txt"],
              ["Larimichthys-crocea_",              "000_OnlyLongestPEP_Larimichthys-crocea_rna.gbk.txt",                                                  "000_OnlyLongestCDS_Larimichthys-crocea_rna.gbk.txt"],
              ["Sparus-aurata_",                    "000_OnlyLongestPEP_Sparus_aurata.fSpaAur1.1.pep.all.fa",                                              "000_OnlyLongestCDS_Sparus_aurata.fSpaAur1.1.cds.all.fa"],
              ["Mola-mola_",                        "000_OnlyLongestPEP_Mola_mola.ASM169857v1.pep.all.fa",                                                 "000_OnlyLongestCDS_Mola_mola.ASM169857v1.cds.all.fa"],
              ["Tetraodon-nigroviridis_",           "000_OnlyLongestPEP_Tetraodon_nigroviridis.TETRAODON8.pep.all.fa",                                     "000_OnlyLongestCDS_Tetraodon_nigroviridis.TETRAODON8.cds.all.fa"],
              ["Takifugu-rubripes_",                "000_OnlyLongestPEP_Takifugu_rubripes.FUGU4.pep.all.fa",                                               "000_OnlyLongestCDS_Takifugu_rubripes.FUGU4.cds.all.fa"],
              ["Takifugu-rubripes-1_",              "000_OnlyLongestPEP_Takifugu-rubripes_GCF_000180615.1_FUGU5_rna.gbff.txt",                             "000_OnlyLongestCDS_Takifugu-rubripes_GCF_000180615.1_FUGU5_rna.gbff.txt"],
              ["Cynoglossus-semilaevis_",           "000_OnlyLongestPEP_Cynoglossus-semilaevis_rna.gbk.txt",                                               "000_OnlyLongestCDS_Cynoglossus-semilaevis_rna.gbk.txt"],
              ["Paralichthys-olivaceus_",           "000_OnlyLongestPEP_Paralichthys-olivaceus_rna.gbk.txt",                                               "000_OnlyLongestCDS_Paralichthys-olivaceus_rna.gbk.txt"],
              ["Scophthalmus-maximus_",             "000_OnlyLongestPEP_Scophthalmus_maximus.ASM318616v1.pep.all.fa",                                      "000_OnlyLongestCDS_Scophthalmus_maximus.ASM318616v1.cds.all.fa"],
              ["Lates-calcarifer_",                 "000_OnlyLongestPEP_Lates-calcarifer_rna.gbk.txt",                                                     "000_OnlyLongestCDS_Lates-calcarifer_rna.gbk.txt"],
              ["Echeneis-naucrates_",               "000_OnlyLongestPEP_Echeneis_naucrates.fEcheNa1.1.pep.all.fa",                                         "000_OnlyLongestCDS_Echeneis_naucrates.fEcheNa1.1.cds.all.fa"],
              ["Seriola-lalandi-dorsalis_",         "000_OnlyLongestPEP_Seriola-lalandi-dorsalis_rna.gbk.txt",                                             "000_OnlyLongestCDS_Seriola-lalandi-dorsalis_rna.gbk.txt"],
              ["Seriola-dumerili_",                 "000_OnlyLongestPEP_Seriola-dumerili_rna.gbk.txt",                                                     "000_OnlyLongestCDS_Seriola-dumerili_rna.gbk.txt"],
              ["Stegastes-partitus_",               "000_OnlyLongestPEP_Stegastes-partitus_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Stegastes-partitus_rna.gbk.txt"],
              ["Amphiprion-ocellaris_",             "000_OnlyLongestPEP_Amphiprion-ocellaris_rna.gbk.txt",                                                 "000_OnlyLongestCDS_Amphiprion-ocellaris_rna.gbk.txt"],
              ["Acanthochromis-polyacanthus_",      "000_OnlyLongestPEP_Acanthochromis-polyacanthus_rna.gbk.txt",                                          "000_OnlyLongestCDS_Acanthochromis-polyacanthus_rna.gbk.txt"],
              ["Gouania-willdenowi_",               "000_OnlyLongestPEP_Gouania_willdenowi.fGouWil2.1.pep.all.fa",                                         "000_OnlyLongestCDS_Gouania_willdenowi.fGouWil2.1.cds.all.fa"],
              ["Salarias-fasciatus_",               "000_OnlyLongestPEP_Salarias_fasciatus.fSalaFa1.1.pep.all.fa",                                         "000_OnlyLongestCDS_Salarias_fasciatus.fSalaFa1.1.cds.all.fa"],
              ["Amphilophus-citrinellus_",          "000_OnlyLongestPEP_Amphilophus_citrinellus.Midas_v5.pep.all.fa",                                      "000_OnlyLongestCDS_Amphilophus_citrinellus.Midas_v5.cds.all.fa"],
              ["Oreochromis-niloticus_",            "000_OnlyLongestPEP_Oreochromis_niloticus.Orenil1.0.pep.all.fa",                                       "000_OnlyLongestCDS_Oreochromis_niloticus.Orenil1.0.cds.all.fa"],
              ["Neolamprologus-brichardi_",         "000_OnlyLongestPEP_Neolamprologus-brichardi_rna.gbk.txt",                                             "000_OnlyLongestCDS_Neolamprologus-brichardi_rna.gbk.txt"],
              ["Pundamilia-nyererei_",              "000_OnlyLongestPEP_Pundamilia-nyererei_rna.gbk.txt",                                                  "000_OnlyLongestCDS_Pundamilia-nyererei_rna.gbk.txt"],
              ["Maylandia-zebra_",                  "000_OnlyLongestPEP_Maylandia-zebra_rna.gbk.txt",                                                      "000_OnlyLongestCDS_Maylandia-zebra_rna.gbk.txt"],
              ["Haplochromis-burtoni_",             "000_OnlyLongestPEP_Haplochromis-burtoni_rna.gbk.txt",                                                 "000_OnlyLongestCDS_Haplochromis-burtoni_rna.gbk.txt"],
              ["Hippocampus-comes_",                "000_OnlyLongestPEP_Hippocampus-comes_rna.gbk.txt",                                                    "000_OnlyLongestCDS_Hippocampus-comes_rna.gbk.txt"],
              ["Sphaeramia-orbicularis_",           "000_OnlyLongestPEP_Sphaeramia_orbicularis.fSphaOr1.1.pep.all.fa",                                     "000_OnlyLongestCDS_Sphaeramia_orbicularis.fSphaOr1.1.cds.all.fa"],
              ["Neogobius-melanostomus_",           "000_OnlyLongestPEP_Neogobius_melanostomus.RGoby_Basel_V2.pep.all.fa",                                 "000_OnlyLongestCDS_Neogobius_melanostomus.RGoby_Basel_V2.cds.all.fa"],
              ["Periophthalmus-magnuspinnatus_",    "000_OnlyLongestPEP_Periophthalmus_magnuspinnatus.PM.fa.pep.all.fa",                                   "000_OnlyLongestCDS_Periophthalmus_magnuspinnatus.PM.fa.cds.all.fa"],
              ["Boleophthalmus-pectinirostris_",    "000_OnlyLongestPEP_Boleophthalmus-pectinirostris_rna.gbk.txt",                                        "000_OnlyLongestCDS_Boleophthalmus-pectinirostris_rna.gbk.txt"],
              ["Monopterus-albus_",                 "000_OnlyLongestPEP_Monopterus-albus_rna.gbk.txt",                                                     "000_OnlyLongestCDS_Monopterus-albus_rna.gbk.txt"],
              ["Mastacembelus-armatus_",            "000_OnlyLongestPEP_Mastacembelus_armatus.fMasArm1.1.pep.all.fa",                                      "000_OnlyLongestCDS_Mastacembelus_armatus.fMasArm1.1.cds.all.fa"],
              ["Anabas-testudineus_",               "000_OnlyLongestPEP_Anabas_testudineus.fAnaTes1.1.pep.all.fa",                                         "000_OnlyLongestCDS_Anabas_testudineus.fAnaTes1.1.cds.all.fa"],
              ["Betta-splendens_",                  "000_OnlyLongestPEP_Betta_splendens.fBetSpl5.2.pep.all.fa",                                            "000_OnlyLongestCDS_Betta_splendens.fBetSpl5.2.cds.all.fa"],
              ["Nothobranchius-furzeri_",           "000_OnlyLongestPEP_Nothobranchius-furzeri_rna.gbk.txt",                                               "000_OnlyLongestCDS_Nothobranchius-furzeri_rna.gbk.txt"],
              ["Austrofundulus-limnaeus_",          "000_OnlyLongestPEP_Austrofundulus-limnaeus_rna.gbk.txt",                                              "000_OnlyLongestCDS_Austrofundulus-limnaeus_rna.gbk.txt"],
              ["Kryptolebias-marmoratus_",          "000_OnlyLongestPEP_Kryptolebias_marmoratus.ASM164957v1.pep.all.fa",                                   "000_OnlyLongestCDS_Kryptolebias_marmoratus.ASM164957v1.cds.all.fa"],
              ["Cyprinodon-variegatus_",            "000_OnlyLongestPEP_Cyprinodon-variegatus_rna.gbk.txt",                                                "000_OnlyLongestCDS_Cyprinodon-variegatus_rna.gbk.txt"],
              ["Fundulus-heteroclitus_",            "000_OnlyLongestPEP_Fundulus-heteroclitus_rna.gbk.txt",                                                "000_OnlyLongestCDS_Fundulus-heteroclitus_rna.gbk.txt"],
              ["Xiphophorus-maculatus_",            "000_OnlyLongestPEP_Xiphophorus_maculatus.Xipmac4.4.2.pep.all.fa",                                     "000_OnlyLongestCDS_Xiphophorus_maculatus.Xipmac4.4.2.cds.all.fa"],
              ["Poecilia-formosa_",                 "000_OnlyLongestPEP_Poecilia_formosa.PoeFor_5.1.2.pep.all.fa",                                         "000_OnlyLongestCDS_Poecilia_formosa.PoeFor_5.1.2.cds.all.fa"],
              ["Poecilia-reticulata_",              "000_OnlyLongestPEP_Poecilia-reticulata_rna.gbk.txt",                                                  "000_OnlyLongestCDS_Poecilia-reticulata_rna.gbk.txt"],
              ["Cololabis-saira_",                  "000_OnlyLongestPEP_Cololabis-sairaGCF_033807715.1_fColSai1.1_rna.gbff.txt",                           "000_OnlyLongestCDS_Cololabis-sairaGCF_033807715.1_fColSai1.1_rna.gbff.txt"],
              ["Oryzias-melastigma_",               "000_OnlyLongestPEP_Oryzias-melastigma_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Oryzias-melastigma_rna.gbk.txt"],
              ["Oryzias-melastigma-E102_",          "000_OnlyLongestPEP_Oryzias_javanicus.OJAV_1.1.pep.all.fa",                                            "000_OnlyLongestCDS_Oryzias_javanicus.OJAV_1.1.cds.all.fa"],
              ["Oryzias-javanicus_",                "000_OnlyLongestPEP_Oryzias_latipes_hni.ASM223471v1.pep.all.fa",                                       "000_OnlyLongestCDS_Oryzias_latipes_hni.ASM223471v1.cds.all.fa"],
              ["Oryzias-sinensis_",                 "000_OnlyLongestPEP_Oryzias_latipes_hsok.ASM223469v1.pep.all.fa",                                      "000_OnlyLongestCDS_Oryzias_latipes_hsok.ASM223469v1.cds.all.fa"],
              ["Oryzias-latipes-HSOK_",             "000_OnlyLongestPEP_Oryzias_latipes.ASM223467v1.pep.all.fa",                                           "000_OnlyLongestCDS_Oryzias_latipes.ASM223467v1.cds.all.fa"],
              ["Oryzias-latipes-HNI_",              "000_OnlyLongestPEP_Oryzias_melastigma.Om_v0.7.RACA.pep.all.fa",                                       "000_OnlyLongestCDS_Oryzias_melastigma.Om_v0.7.RACA.cds.all.fa"],
              ["Oryzias-latipes-E102_",             "000_OnlyLongestPEP_Oryzias_sinensis.ASM858656v1.pep.all.fa",                                          "000_OnlyLongestCDS_Oryzias_sinensis.ASM858656v1.cds.all.fa"],

              ["Oryzias-latipes_",                  "000_OnlyLongestPEP_Oryzias_latipes.MEDAKA1.pep.all.fa",                                               "000_OnlyLongestCDS_Oryzias_latipes.MEDAKA1.cds.all.fa"],
              ["Oryzias-latipes-1_",                "000_OnlyLongestPEP_Oryzias-latipes_GCF_002234675.1_ASM223467v1_rna.gbff.txt",                         "000_OnlyLongestCDS_Oryzias-latipes_GCF_002234675.1_ASM223467v1_rna.gbff.txt"],

              ### Sarcopterygii
              ["Latimeria-chalumnae_",              "000_OnlyLongestPEP_Latimeria_chalumnae.LatCha1.pep.all.fa",                                           "000_OnlyLongestCDS_Latimeria_chalumnae.LatCha1.cds.all.fa"],
              ["Latimeria-chalumnae-1_",            "000_OnlyLongestPEP_Latimeria-chalumnae_GCF_000225785.1_LatCha1_rna.gbff.txt",                         "000_OnlyLongestCDS_Latimeria-chalumnae_GCF_000225785.1_LatCha1_rna.gbff.txt"],
              
              ### Amphibia
              #["Rhinatrema-bivittatum_",            "005_Rhinatrema-bivittatum_longest_orfs.pep",                                                          "005_Rhinatrema-bivittatum_longest_orfs.cds"],
              ["Rhinatrema-bivittatum_",            "000_OnlyLongestPEP_prot_Rhinatrema-bivittatum_GCF_901001135.1_aRhiBiv1.1_rna.gbff.txt",               "000_OnlyLongestCDS_nucl_Rhinatrema-bivittatum_GCF_901001135.1_aRhiBiv1.1_rna.gbff.txt"],
              ["Caecilia-tentaculata_",             "005_Caecilia-tentaculata_longest_orfs.pep",                                                           "005_Caecilia-tentaculata_longest_orfs.cds"],
              ["Typhlonectes-compressicauda_",      "005_Typhlonectes-compressicauda_longest_orfs.pep",                                                    "005_Typhlonectes-compressicauda_longest_orfs.cds"],
              ["Microcaecilia-dermatophaga_",       "005_Microcaecilia-dermatophaga_longest_orfs.pep",                                                     "005_Microcaecilia-dermatophaga_longest_orfs.cds"],
              ["Microcaecilia-unicolor_",           "000_OnlyLongestPEP_prot_Microcaecilia-unicolor_GCF_901765095.1_aMicUni1.1_rna.gbff.txt",              "000_OnlyLongestCDS_nucl_Microcaecilia-unicolor_GCF_901765095.1_aMicUni1.1_rna.gbff.txt"],
              #["Microcaecilia-unicolor_",           "005_Microcaecilia-unicolor_longest_orfs.pep",                                                         "005_Microcaecilia-unicolor_longest_orfs.cds"],
              ["Tylototriton-wenxianensis_",        "005_Tylototriton-wenxianensis_longest_orfs.pep",                                                      "005_Tylototriton-wenxianensis_longest_orfs.cds"],
              ["Nanorana-parkeri_",                 "000_OnlyLongestPEP_Nanorana-parkeri_GCF_000935625.1_ASM93562v1_rna.gbff.txt",                         "000_OnlyLongestCDS_Nanorana-parkeri_GCF_000935625.1_ASM93562v1_rna.gbff.txt"],
              ["Xenopus-tropicalis_",               "000_OnlyLongestPEP_Xenopus_tropicalis.JGI_4.2.pep.all.fa",                                            "000_OnlyLongestCDS_Xenopus_tropicalis.JGI_4.2.cds.all.fa"],
              ["Xenopus-tropicalis-1_",             "000_OnlyLongestPEP_Xenopus-tropicalis_frog.1.rna.gbff.txt",                                           "000_OnlyLongestCDS_Xenopus-tropicalis_frog.1.rna.gbff.txt"],
              ["Xenopus-laevis_",                   "000_OnlyLongestPEP_Xenopus-laevis_rna.gbk.txt",                                                       "000_OnlyLongestCDS_Xenopus-laevis_rna.gbk.txt"],
              ["Xenopus-laevis-1_",                 "000_OnlyLongestPEP_GCF_001663975.1_Xenopus_laevis_v2_rna.gbff.txt",                                   "000_OnlyLongestCDS_GCF_001663975.1_Xenopus_laevis_v2_rna.gbff.txt"],
              
              #### Lepidosauria
              ["Sphenodon-punctatus_",              "000_OnlyLongestPEP_Sphenodon_punctatus.ASM311381v1.pep.all.fa",                                       "000_OnlyLongestCDS_Sphenodon_punctatus.ASM311381v1.cds.all.fa"],
              ["Gekko-japonicus_",                  "000_OnlyLongestPEP_Gekko-japonicus_rna.gbk.txt",                                                      "000_OnlyLongestCDS_Gekko-japonicus_rna.gbk.txt"],
              ["Podarcis-muralis_",                 "000_OnlyLongestPEP_prot_Podarcis-muralis_GCF_004329235.1_PodMur_1.0_rna.gbff.txt",                    "000_OnlyLongestCDS_nucl_Podarcis-muralis_GCF_004329235.1_PodMur_1.0_rna.gbff.txt"],
              ["Salvator-merianae_",                "000_OnlyLongestPEP_Salvator_merianae.HLtupMer3.pep.all.fa",                                           "000_OnlyLongestCDS_Salvator_merianae.HLtupMer3.cds.all.fa"],
              ["Pogona-vitticeps_",                 "000_OnlyLongestPEP_Pogona-vitticeps_rna.gbk.txt",                                                     "000_OnlyLongestCDS_Pogona-vitticeps_rna.gbk.txt"],
              ["Anolis-carolinensis_",              "000_OnlyLongestPEP_Anolis_carolinensis.AnoCar2.0.pep.all.fa",                                         "000_OnlyLongestCDS_Anolis_carolinensis.AnoCar2.0.cds.all.fa"],
              ["Anolis-carolinensis-1_",            "000_OnlyLongestPEP_Anolis-carolinensis_GCF_000090745.1_AnoCar2.0_rna.gbff.txt",                       "000_OnlyLongestCDS_Anolis-carolinensis_GCF_000090745.1_AnoCar2.0_rna.gbff.txt"],
              ["Python-bivittatus_",                "000_OnlyLongestPEP_Python-bivittatus_rna.gbk.txt",                                                    "000_OnlyLongestCDS_Python-bivittatus_rna.gbk.txt"],
              ["Protobothrops-mucrosquamatus_",     "000_OnlyLongestPEP_Protobothrops-mucrosquamatus_rna.gbk.txt",                                         "000_OnlyLongestCDS_Protobothrops-mucrosquamatus_rna.gbk.txt"],
              ["Protobothrops-flavoviridis_",       "000_OnlyLongestPEP_Protobothrops-flavoviridis_habu1_augustus.aa",                                     "000_OnlyLongestCDS_Protobothrops-flavoviridis_habu1_augustus.fa"],
              ["Pantherophis-guttatus_",            "000_OnlyLongestPEP_prot_GCF_001185365.1_UNIGE_PanGut_3.0_rna.gbff.txt",                               "000_OnlyLongestCDS_nucl_GCF_001185365.1_UNIGE_PanGut_3.0_rna.gbff.txt"],
              ["Thamnophis-sirtalis_",              "000_OnlyLongestPEP_Thamnophis-sirtalis_rna.gbk.txt",                                                  "000_OnlyLongestCDS_Thamnophis-sirtalis_rna.gbk.txt"],
              ["Thamnophis-elegans_",               "000_OnlyLongestPEP_prot_GCF_009769535.1_rThaEle1.pri_rna.gbff.txt",                                   "000_OnlyLongestCDS_nucl_GCF_009769535.1_rThaEle1.pri_rna.gbff.txt"],
              ["Naja-naja_",                        "000_OnlyLongestPEP_Naja_naja.Nana_v5.pep.all.fa",                                                     "000_OnlyLongestCDS_Naja_naja.Nana_v5.cds.all.fa"],
              ["Laticauda-laticaudata_",            "000_OnlyLongestPEP_Laticauda_laticaudata.latLat_1.0.pep.all.fa",                                      "000_OnlyLongestCDS_Laticauda_laticaudata.latLat_1.0.cds.all.fa"],
              ["Pseudonaja-textilis_",              "000_OnlyLongestPEP_prot_Pseudonaja-textilis_GCF_900518735.1_EBS10Xv2-PRI_rna.gbff.txt",               "000_OnlyLongestCDS_nucl_Pseudonaja-textilis_GCF_900518735.1_EBS10Xv2-PRI_rna.gbff.txt"],
              ["Notechis-scutatus_",                "000_OnlyLongestPEP_prot_Notechis-scutatus_GCF_900518725.1_TS10Xv2-PRI_rna.gbff.txt",                  "000_OnlyLongestCDS_nucl_Notechis-scutatus_GCF_900518725.1_TS10Xv2-PRI_rna.gbff.txt"],
              
              ### Archosauria
              ["Pelodiscus-sinensis_",              "000_OnlyLongestPEP_Pelodiscus_sinensis.PelSin_1.0.pep.all.fa",                                        "000_OnlyLongestCDS_Pelodiscus_sinensis.PelSin_1.0.cds.all.fa"],
              ["Chelonia-mydas_",                   "000_OnlyLongestPEP_Chelonia-mydas_rna.gbk.txt",                                                       "000_OnlyLongestCDS_Chelonia-mydas_rna.gbk.txt"],
              ["Chelonoidis-abingdonii_",           "000_OnlyLongestPEP_Chelonoidis_abingdonii.ASM359739v1.pep.all.fa",                                    "000_OnlyLongestCDS_Chelonoidis_abingdonii.ASM359739v1.cds.all.fa"],
              ["Terrapene-mexicana_",               "000_OnlyLongestPEP_Terrapene-mexicana_GCF_002925995.1_T_m_triunguis-1.0_rna.gbff.txt",                "000_OnlyLongestCDS_Terrapene-mexicana_GCF_002925995.1_T_m_triunguis-1.0_rna.gbff.txt"],
              ["Chrysemys-picta_",                  "000_OnlyLongestPEP_Chrysemys-picta_rna.gbk.txt",                                                      "000_OnlyLongestCDS_Chrysemys-picta_rna.gbk.txt"],
              ["Crocodylus-porosus_",               "000_OnlyLongestPEP_Crocodylus-porosus_GCF_001723895.1_CroPor_comp1_rna.gbff.txt",                     "000_OnlyLongestCDS_Crocodylus-porosus_GCF_001723895.1_CroPor_comp1_rna.gbff.txt"],
              ["Gavialis-gangeticus_",              "000_OnlyLongestPEP_Gavialis-gangeticus_GCF_001723915.1_GavGan_comp1_rna.gbff.txt",                    "000_OnlyLongestCDS_Gavialis-gangeticus_GCF_001723915.1_GavGan_comp1_rna.gbff.txt"],
              ["Alligator-mississippiensis_",       "000_OnlyLongestPEP_Alligator-mississippiensis_rna.gbk.txt",                                           "000_OnlyLongestCDS_Alligator-mississippiensis_rna.gbk.txt"],
              ["Alligator-sinensis_",               "000_OnlyLongestPEP_Alligator-sinensis_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Alligator-sinensis_rna.gbk.txt"],

              ## Aves
              ["Struthio-camelus_",                 "000_OnlyLongestPEP_Struthio-camelus_GCF_000698965.1_ASM69896v1_rna.gbff.txt",                         "000_OnlyLongestCDS_Struthio-camelus_GCF_000698965.1_ASM69896v1_rna.gbff.txt"],
              ["Tinamus-guttatus_",                 "000_OnlyLongestPEP_Tinamus-guttatus_GCF_000705375.1_ASM70537v2_rna.gbff.txt",                         "000_OnlyLongestCDS_Tinamus-guttatus_GCF_000705375.1_ASM70537v2_rna.gbff.txt"],
              ["Nothoprocta-perdicaria_",           "000_OnlyLongestPEP_Nothoprocta-perdicaria_GCF_003342845.1_notPer1_rna.gbff.txt",                      "000_OnlyLongestCDS_Nothoprocta-perdicaria_GCF_003342845.1_notPer1_rna.gbff.txt"],
              ["Apteryx-rowi_",                     "000_OnlyLongestPEP_Apteryx-rowi_GCF_003343035.1_aptRow1_rna.gbff.txt",                                "000_OnlyLongestCDS_Apteryx-rowi_GCF_003343035.1_aptRow1_rna.gbff.txt"],
              ["Apteryx-australis_",                "000_OnlyLongestPEP_Apteryx-australis_GCF_001039765.1_AptMant0_rna.gbff.txt",                          "000_OnlyLongestCDS_Apteryx-australis_GCF_001039765.1_AptMant0_rna.gbff.txt"],
              ["Dromaius-novaehollandiae_",         "000_OnlyLongestPEP_Dromaius-novaehollandiae_GCF_003342905.1_droNov1_rna.gbff.txt",                    "000_OnlyLongestCDS_Dromaius-novaehollandiae_GCF_003342905.1_droNov1_rna.gbff.txt"],

              ["Anser-cygnoides_",                  "000_OnlyLongestPEP_Anser-cygnoides_GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_rna.gbff.txt",             "000_OnlyLongestCDS_Anser-cygnoides_GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_rna.gbff.txt"],
              ["Anas-platyrhynchos_",               "000_OnlyLongestPEP_Anas_platyrhynchos.BGI_duck_1.0.pep.all.fa",                                       "000_OnlyLongestCDS_Anas_platyrhynchos.BGI_duck_1.0.cds.all.fa"],
              ["Anas-platyrhynchos-1_",             "000_OnlyLongestPEP_Anas-platyrhynchos_GCF_000355885.1_BGI_duck_1.0_rna.gbff.txt",                     "000_OnlyLongestCDS_Anas-platyrhynchos_GCF_000355885.1_BGI_duck_1.0_rna.gbff.txt"],
              ["Meleagris-gallopavo_",              "000_OnlyLongestPEP_Meleagris_gallopavo.UMD2.pep.all.fa",                                              "000_OnlyLongestCDS_Meleagris_gallopavo.UMD2.cds.all.fa"],
              ["Meleagris-gallopavo-1_",            "000_OnlyLongestPEP_Meleagris-gallopavo_GCF_000146605.2_Turkey_5.0_rna.gbff.txt",                      "000_OnlyLongestCDS_Meleagris-gallopavo_GCF_000146605.2_Turkey_5.0_rna.gbff.txt"],
              ["Coturnix-japonica_",                "000_OnlyLongestPEP_Coturnix-japonica_GCF_001577835.1_Coturnix_japonica_2.0_rna.gbff.txt",             "000_OnlyLongestCDS_Coturnix-japonica_GCF_001577835.1_Coturnix_japonica_2.0_rna.gbff.txt"],
              ["Gallus-gallus-E104_",               "000_OnlyLongestPEP_Gallus_gallus-E104.GRCg6a.pep.all.fa",                                             "000_OnlyLongestCDS_Gallus_gallus-E104.GRCg6a.cds.all.fa"],
              #["Gallus-gallus-E102_",               "000_OnlyLongestPEP_Gallus_gallus-E102.GRCg6a.pep.all.fa",                                             "000_OnlyLongestCDS_Gallus_gallus-E102.GRCg6a.cds.all.fa"],
              #["Gallus-gallus-E99_",                "000_OnlyLongestPEP_Gallus_gallus.GRCg6a.pep.all.fa",                                                  "000_OnlyLongestCDS_Gallus_gallus.GRCg6a.cds.all.fa"],
              #["Gallus-gallus_",                    "000_OnlyLongestPEP_Gallus_gallus.Gallus_gallus-5.0.pep.all.fa",                                       "000_OnlyLongestCDS_Gallus_gallus.Gallus_gallus-5.0.cds.all.fa"],
              #["Gallus-gallus-1_",                  "000_OnlyLongestPEP_Gallus-gallus_GCF_000002315.5_GRCg6a_rna.gbff.txt",                                "000_OnlyLongestCDS_Gallus-gallus_GCF_000002315.5_GRCg6a_rna.gbff.txt"],
              ["Gallus-gallus-RS2_",                "000_OnlyLongestPEP_Gallus-gallus-RS2_GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_rna.gbff.txt",       "000_OnlyLongestCDS_Gallus-gallus-RS2_GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_rna.gbff.txt"],
              ["Columba-livia_",                    "000_OnlyLongestPEP_Columba-livia_GCF_000337935.1_Cliv_1.0_rna.gbff.txt",                              "000_OnlyLongestCDS_Columba-livia_GCF_000337935.1_Cliv_1.0_rna.gbff.txt"],
              ["Pterocles-gutturalis_",             "000_OnlyLongestPEP_Pterocles_gutturalis_GCF_000699245.1_ASM69924v1_rna.gbff.txt",                     "000_OnlyLongestCDS_Pterocles_gutturalis_GCF_000699245.1_ASM69924v1_rna.gbff.txt"],
              ["Tauraco-erythrolophus_",            "000_OnlyLongestPEP_Tauraco_erythrolophus_GCF_000709365.1_ASM70936v1_rna.gbff.txt",                    "000_OnlyLongestCDS_Tauraco_erythrolophus_GCF_000709365.1_ASM70936v1_rna.gbff.txt"],
              ["Cuculus-canorus_",                  "000_OnlyLongestPEP_Cuculus_canorus_GCF_000709325.1_ASM70932v1_rna.gbff.txt",                          "000_OnlyLongestCDS_Cuculus_canorus_GCF_000709325.1_ASM70932v1_rna.gbff.txt"],
              ["Antrostomus-carolinensis_",         "000_OnlyLongestPEP_Antrostomus_carolinensis_GCF_000700745.1_ASM70074v1_rna.gbff.txt",                 "000_OnlyLongestCDS_Antrostomus_carolinensis_GCF_000700745.1_ASM70074v1_rna.gbff.txt"],
              ["Calypte-anna_",                     "000_OnlyLongestPEP_Calypte_anna_GCF_000699085.1_ASM69908v1_rna.gbff.txt",                             "000_OnlyLongestCDS_Calypte_anna_GCF_000699085.1_ASM69908v1_rna.gbff.txt"],
              ["Chaetura-pelagica_",                "000_OnlyLongestPEP_Chaetura_pelagica_GCF_000747805.1_ChaPel_1.0_rna.gbff.txt",                        "000_OnlyLongestCDS_Chaetura_pelagica_GCF_000747805.1_ChaPel_1.0_rna.gbff.txt"],
              ["Opisthocomus-hoazin_",              "000_OnlyLongestPEP_Opisthocomus_hoazin_GCF_000692075.1_ASM69207v1_rna.gbff.txt",                      "000_OnlyLongestCDS_Opisthocomus_hoazin_GCF_000692075.1_ASM69207v1_rna.gbff.txt"],
              ["Chlamydotis-macqueenii_",           "000_OnlyLongestPEP_Chlamydotis_macqueenii_GCF_000695195.1_ASM69519v1_rna.gbff.txt",                   "000_OnlyLongestCDS_Chlamydotis_macqueenii_GCF_000695195.1_ASM69519v1_rna.gbff.txt"],
              ["Mesitornis-unicolor_",              "000_OnlyLongestPEP_Mesitornis_unicolor_GCF_000695765.1_ASM69576v1_rna.gbff.txt",                      "000_OnlyLongestCDS_Mesitornis_unicolor_GCF_000695765.1_ASM69576v1_rna.gbff.txt"],
              ["Eurypyga-helias_",                  "000_OnlyLongestPEP_Eurypyga_helias_GCF_000690775.1_ASM69077v1_rna.gbff.txt",                          "000_OnlyLongestCDS_Eurypyga_helias_GCF_000690775.1_ASM69077v1_rna.gbff.txt"],
              ["Balearica-regulorum_",              "000_OnlyLongestPEP_Balearica_regulorum_GCF_000709895.1_ASM70989v1_rna.gbff.txt",                      "000_OnlyLongestCDS_Balearica_regulorum_GCF_000709895.1_ASM70989v1_rna.gbff.txt"],
              ["Calidris-pugnax_",                  "000_OnlyLongestPEP_Calidris_pugnax_GCF_001431845.1_ASM143184v1_rna.gbff.txt",                         "000_OnlyLongestCDS_Calidris_pugnax_GCF_001431845.1_ASM143184v1_rna.gbff.txt"],
              ["Charadrius-vociferus_",             "000_OnlyLongestPEP_Charadrius_vociferus_GCF_000708025.1_ASM70802v2_rna.gbff.txt",                     "000_OnlyLongestCDS_Charadrius_vociferus_GCF_000708025.1_ASM70802v2_rna.gbff.txt"],
              ["Gavia-stellata_",                   "000_OnlyLongestPEP_Gavia_stellata_GCF_000690875.1_ASM69087v1_rna.gbff.txt",                           "000_OnlyLongestCDS_Gavia_stellata_GCF_000690875.1_ASM69087v1_rna.gbff.txt"],
              ["Pygoscelis-adeliae_",               "000_OnlyLongestPEP_Pygoscelis_adeliae_GCF_000699105.1_ASM69910v1_rna.gbff.txt",                       "000_OnlyLongestCDS_Pygoscelis_adeliae_GCF_000699105.1_ASM69910v1_rna.gbff.txt"],
              ["Aptenodytes-forsteri_",             "000_OnlyLongestPEP_Aptenodytes_forsteri_GCF_000699145.1_ASM69914v1_rna.gbff.txt",                     "000_OnlyLongestCDS_Aptenodytes_forsteri_GCF_000699145.1_ASM69914v1_rna.gbff.txt"],
              ["Fulmarus-glacialis_",               "000_OnlyLongestPEP_Fulmarus_glacialis_GCF_000690835.1_ASM69083v1_rna.gbff.txt",                       "000_OnlyLongestCDS_Fulmarus_glacialis_GCF_000690835.1_ASM69083v1_rna.gbff.txt"],
              ["Phaethon-lepturus_",                "000_OnlyLongestPEP_Phaethon_lepturus_GCF_000687285.1_ASM68728v1_rna.gbff.txt",                        "000_OnlyLongestCDS_Phaethon_lepturus_GCF_000687285.1_ASM68728v1_rna.gbff.txt"],
              ["Phalacrocorax-carbo_",              "000_OnlyLongestPEP_Phalacrocorax_carbo_GCF_000708925.1_ASM70892v1_rna.gbff.txt",                      "000_OnlyLongestCDS_Phalacrocorax_carbo_GCF_000708925.1_ASM70892v1_rna.gbff.txt"],
              ["Pelecanus-crispus_",                "000_OnlyLongestPEP_Pelecanus_crispus_GCF_000687375.1_ASM68737v1_rna.gbff.txt",                        "000_OnlyLongestCDS_Pelecanus_crispus_GCF_000687375.1_ASM68737v1_rna.gbff.txt"],
              ["Nipponia-nippon_",                  "000_OnlyLongestPEP_Nipponia_nippon_GCF_000708225.1_ASM70822v1_rna.gbff.txt",                          "000_OnlyLongestCDS_Nipponia_nippon_GCF_000708225.1_ASM70822v1_rna.gbff.txt"],
              ["Egretta-garzetta_",                 "000_OnlyLongestPEP_Egretta_garzetta_GCF_000687185.1_ASM68718v1_rna.gbff.txt",                         "000_OnlyLongestCDS_Egretta_garzetta_GCF_000687185.1_ASM68718v1_rna.gbff.txt"],
              ["Athene-cunicularia_",               "000_OnlyLongestPEP_Athene_cunicularia_GCF_003259725.1_athCun1_rna.gbff.txt",                          "000_OnlyLongestCDS_Athene_cunicularia_GCF_003259725.1_athCun1_rna.gbff.txt"],
              ["Tyto-alba_",                        "000_OnlyLongestPEP_Tyto_alba_GCF_000687205.1_ASM68720v1_rna.gbff.txt",                                "000_OnlyLongestCDS_Tyto_alba_GCF_000687205.1_ASM68720v1_rna.gbff.txt"],
              ["Colius-striatus_",                  "000_OnlyLongestPEP_Colius_striatus_GCF_000690715.1_ASM69071v1_rna.gbff.txt",                          "000_OnlyLongestCDS_Colius_striatus_GCF_000690715.1_ASM69071v1_rna.gbff.txt"],
              ["Apaloderma-vittatum_",              "000_OnlyLongestPEP_Apaloderma_vittatum_GCF_000703405.1_ASM70340v1_rna.gbff.txt",                      "000_OnlyLongestCDS_Apaloderma_vittatum_GCF_000703405.1_ASM70340v1_rna.gbff.txt"],
              ["Buceros-rhinoceros_",               "000_OnlyLongestPEP_Buceros_rhinoceros_GCF_000710305.1_ASM71030v1_rna.gbff.txt",                       "000_OnlyLongestCDS_Buceros_rhinoceros_GCF_000710305.1_ASM71030v1_rna.gbff.txt"],
              ["Leptosomus-discolor_",              "000_OnlyLongestPEP_Leptosomus_discolor_GCF_000691785.1_ASM69178v1_rna.gbff.txt",                      "000_OnlyLongestCDS_Leptosomus_discolor_GCF_000691785.1_ASM69178v1_rna.gbff.txt"],
              ["Merops-nubicus_",                   "000_OnlyLongestPEP_Merops_nubicus_GCF_000691845.1_ASM69184v1_rna.gbff.txt",                           "000_OnlyLongestCDS_Merops_nubicus_GCF_000691845.1_ASM69184v1_rna.gbff.txt"],
              ["Picoides-pubescens_",               "000_OnlyLongestPEP_Picoides_pubescens_GCF_000699005.1_ASM69900v1_rna.gbff.txt",                       "000_OnlyLongestCDS_Picoides_pubescens_GCF_000699005.1_ASM69900v1_rna.gbff.txt"],
              ["Cariama-cristata_",                 "000_OnlyLongestPEP_Cariama_cristata_GCF_000690535.1_ASM69053v1_rna.gbff.txt",                         "000_OnlyLongestCDS_Cariama_cristata_GCF_000690535.1_ASM69053v1_rna.gbff.txt"],
              ["Falco-cherrug_",                    "000_OnlyLongestPEP_Falco_cherrug_GCF_000337975.1_F_cherrug_v1.0_rna.gbff.txt",                        "000_OnlyLongestCDS_Falco_cherrug_GCF_000337975.1_F_cherrug_v1.0_rna.gbff.txt"],
              ["Falco-peregrinus_",                 "000_OnlyLongestPEP_Falco_peregrinus_GCF_000337955.1_F_peregrinus_v1.0_rna.gbff.txt",                  "000_OnlyLongestCDS_Falco_peregrinus_GCF_000337955.1_F_peregrinus_v1.0_rna.gbff.txt"],
              ["Aquila-chrysaetos_",                "000_OnlyLongestPEP_GCF_000766835.1_Aquila_chrysaetos-1.0.2_rna.gbff.txt",                             "000_OnlyLongestCDS_GCF_000766835.1_Aquila_chrysaetos-1.0.2_rna.gbff.txt"],
              ["Haliaeetus-leucocephalus_",         "000_OnlyLongestPEP_GCF_000737465.1_Haliaeetus_leucocephalus-4.0_rna.gbff.txt",                        "000_OnlyLongestCDS_GCF_000737465.1_Haliaeetus_leucocephalus-4.0_rna.gbff.txt"],
              ["Haliaeetus-albicilla_",             "000_OnlyLongestPEP_Haliaeetus_albicilla_GCF_000691405.1_ASM69140v1_rna.gbff.txt",                     "000_OnlyLongestCDS_Haliaeetus_albicilla_GCF_000691405.1_ASM69140v1_rna.gbff.txt"],
              ["Nestor-notabilis_",                 "000_OnlyLongestPEP_Nestor_notabilis_GCF_000696875.1_ASM69687v1_rna.gbff.txt",                         "000_OnlyLongestCDS_Nestor_notabilis_GCF_000696875.1_ASM69687v1_rna.gbff.txt"],
              ["Melopsittacus-undulatus_",          "000_OnlyLongestPEP_GCF_000238935.1_Melopsittacus_undulatus_6.3_rna.gbff.txt",                         "000_OnlyLongestCDS_GCF_000238935.1_Melopsittacus_undulatus_6.3_rna.gbff.txt"],
              ["Acanthisitta-chloris_",             "000_OnlyLongestPEP_Acanthisitta_chloris_GCF_000695815.1_ASM69581v1_rna.gbff.txt",                     "000_OnlyLongestCDS_Acanthisitta_chloris_GCF_000695815.1_ASM69581v1_rna.gbff.txt"],
              ["Corvus-cornix_",                    "000_OnlyLongestPEP_Corvus_cornix_GCF_000738735.2_ASM73873v2_rna.gbff.txt",                            "000_OnlyLongestCDS_Corvus_cornix_GCF_000738735.2_ASM73873v2_rna.gbff.txt"],
              ["Corvus-brachyrhynchos_",            "000_OnlyLongestPEP_Corvus_brachyrhynchos_GCF_000691975.1_ASM69197v1_rna.gbff.txt",                    "000_OnlyLongestCDS_Corvus_brachyrhynchos_GCF_000691975.1_ASM69197v1_rna.gbff.txt"],
              ["Cyanistes-caeruleus_",              "000_OnlyLongestPEP_Cyanistes_caeruleus_GCF_002901205.1_cyaCae2_rna.gbff.txt",                         "000_OnlyLongestCDS_Cyanistes_caeruleus_GCF_002901205.1_cyaCae2_rna.gbff.txt"],
              ["Pseudopodoces-humilis_",            "000_OnlyLongestPEP_Pseudopodoces_humilis_GCF_000331425.1_PseHum1.0_rna.gbff.txt",                     "000_OnlyLongestCDS_Pseudopodoces_humilis_GCF_000331425.1_PseHum1.0_rna.gbff.txt"],
              ["Parus-major_",                      "000_OnlyLongestPEP_Parus_major_GCF_001522545.2_Parus_major1.1_rna.gbff.txt",                          "000_OnlyLongestCDS_Parus_major_GCF_001522545.2_Parus_major1.1_rna.gbff.txt"],
              ["Ficedula-albicollis_",              "000_OnlyLongestPEP_Ficedula_albicollis.FicAlb_1.4.pep.all.fa",                                        "000_OnlyLongestCDS_Ficedula_albicollis.FicAlb_1.4.cds.all.fa"],
              ["Ficedula-albicollis-1_",            "000_OnlyLongestPEP_Ficedula_albicollis_GCF_000247815.1_FicAlb1.5_rna.gbff.txt",                       "000_OnlyLongestCDS_Ficedula_albicollis_GCF_000247815.1_FicAlb1.5_rna.gbff.txt"],
              ["Sturnus-vulgaris_",                 "000_OnlyLongestPEP_GCF_001447265.1_Sturnus_vulgaris-1.0_rna.gbff.txt",                                "000_OnlyLongestCDS_GCF_001447265.1_Sturnus_vulgaris-1.0_rna.gbff.txt"],
              ["Serinus-canaria_",                  "000_OnlyLongestPEP_Serinus_canaria_GCF_000534875.1_SCA1_rna.gbff.txt",                                "000_OnlyLongestCDS_Serinus_canaria_GCF_000534875.1_SCA1_rna.gbff.txt"],
              ["Manacus-vitellinus_",               "000_OnlyLongestPEP_Manacus_vitellinus_GCF_001715985.2_ASM171598v2_rna.gbff.txt",                      "000_OnlyLongestCDS_Manacus_vitellinus_GCF_001715985.2_ASM171598v2_rna.gbff.txt"],
              ["Lepidothrix-coronata_",             "000_OnlyLongestPEP_GCF_001604755.1_Lepidothrix_coronata-1.0_rna.gbff.txt",                            "000_OnlyLongestCDS_GCF_001604755.1_Lepidothrix_coronata-1.0_rna.gbff.txt"],
              ["Zonotrichia-albicollis_",           "000_OnlyLongestPEP_GCF_000385455.1_Zonotrichia_albicollis-1.0.1_rna.gbff.txt",                        "000_OnlyLongestCDS_GCF_000385455.1_Zonotrichia_albicollis-1.0.1_rna.gbff.txt"],
              ["Geospiza-fortis_",                  "000_OnlyLongestPEP_Geospiza_fortis_GCF_000277835.1_GeoFor_1.0_rna.gbff.txt",                          "000_OnlyLongestCDS_Geospiza_fortis_GCF_000277835.1_GeoFor_1.0_rna.gbff.txt"],
              ["Taeniopygia-guttata_",              "000_OnlyLongestPEP_Taeniopygia_guttata.taeGut3.2.4.pep.all.fa",                                       "000_OnlyLongestCDS_Taeniopygia_guttata.taeGut3.2.4.cds.all.fa"],
              ["Taeniopygia-guttata-1_",            "000_OnlyLongestPEP_GCF_000151805.1_Taeniopygia_guttata-3.2.4_rna.gbff.txt",                           "000_OnlyLongestCDS_GCF_000151805.1_Taeniopygia_guttata-3.2.4_rna.gbff.txt"],
              ["Lonchura-striata_",                 "000_OnlyLongestPEP_Lonchura_striata_GCF_002197715.1_LonStrDom1_rna.gbff.txt",                         "000_OnlyLongestCDS_Lonchura_striata_GCF_002197715.1_LonStrDom1_rna.gbff.txt"],

              ## Mammalia
              
              ["Ornithorhynchus-anatinus_",         "000_OnlyLongestPEP_Ornithorhynchus_anatinus.OANA5.pep.all.fa",                                        "000_OnlyLongestCDS_Ornithorhynchus_anatinus.OANA5.cds.all.fa"],
              
              ## Marsupialia
              ["Monodelphis-domestica_",            "000_OnlyLongestPEP_Monodelphis_domestica.monDom5.pep.all.fa",                                         "000_OnlyLongestCDS_Monodelphis_domestica.monDom5.cds.all.fa"],
              ["Sarcophilus-harrisii_",             "000_OnlyLongestPEP_Sarcophilus_harrisii.DEVIL7.0.pep.all.fa",                                         "000_OnlyLongestCDS_Sarcophilus_harrisii.DEVIL7.0.cds.all.fa"],
              ["Notamacropus-eugenii_",             "000_OnlyLongestPEP_Notamacropus_eugenii.Meug_1.0.pep.all.fa",                                         "000_OnlyLongestCDS_Notamacropus_eugenii.Meug_1.0.cds.all.fa"],
              ["Vombatus-ursinus_",                 "000_OnlyLongestPEP_Vombatus_ursinus.bare-nosed_wombat_genome_assembly.pep.all.fa",                    "000_OnlyLongestCDS_Vombatus_ursinus.bare-nosed_wombat_genome_assembly.cds.all.fa"],
              ["Phascolarctos-cinereus_",           "000_OnlyLongestPEP_Phascolarctos_cinereus.phaCin_unsw_v4.1.pep.all.fa",                               "000_OnlyLongestCDS_Phascolarctos_cinereus.phaCin_unsw_v4.1.cds.all.fa"],

              ### Xenarthra
              ["Dasypus-novemcinctus_",             "000_OnlyLongestPEP_Dasypus_novemcinctus.Dasnov3.0.pep.all.fa",                                        "000_OnlyLongestCDS_Dasypus_novemcinctus.Dasnov3.0.cds.all.fa"],
              ["Choloepus-hoffmanni_",              "000_OnlyLongestPEP_Choloepus_hoffmanni.choHof1.pep.all.fa",                                           "000_OnlyLongestCDS_Choloepus_hoffmanni.choHof1.cds.all.fa"],
              
              ### Afrotheria
              ["Elephantulus-edwardii_",            "000_OnlyLongestPEP_prot_Elephantulus-edwardii_GCF_000299155.1_EleEdw1.0_rna.gbff.txt",                "000_OnlyLongestCDS_nucl_Elephantulus-edwardii_GCF_000299155.1_EleEdw1.0_rna.gbff.txt"], 
              ["Orycteropus-afer_",                 "000_OnlyLongestPEP_Orycteropus-afer_rna.gbk.txt",                                                     "000_OnlyLongestCDS_Orycteropus-afer_rna.gbk.txt"],
              ["Chrysochloris-asiatica_",           "000_OnlyLongestPEP_prot_Chrysochloris-asiatica_GCF_000296735.1_ChrAsi1.0_rna.gbff.txt",               "000_OnlyLongestCDS_nucl_Chrysochloris-asiatica_GCF_000296735.1_ChrAsi1.0_rna.gbff.txt"], 
              ["Echinops-telfairi_",                "000_OnlyLongestPEP_Echinops-telfairi.TENREC.pep.all.fa",                                              "000_OnlyLongestCDS_Echinops-telfairi.TENREC.cds.all.fa"],
              ["Procavia-capensis_",                "000_OnlyLongestPEP_Procavia_capensis.proCap1.pep.all.fa",                                             "000_OnlyLongestCDS_Procavia_capensis.proCap1.cds.all.fa"],
              ["Loxodonta-africana_",               "000_OnlyLongestPEP_Loxodonta_africana.loxAfr3.pep.all.fa",                                            "000_OnlyLongestCDS_Loxodonta_africana.loxAfr3.cds.all.fa"],
              ["Trichechus-manatus_",               "000_OnlyLongestPEP_Trichechus-manatus_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Trichechus-manatus_rna.gbk.txt"],
            
              ####Laurasiatheria            
              ["Erinaceus-europaeus_",              "000_OnlyLongestPEP_Erinaceus_europaeus.HEDGEHOG.pep.all.fa",                                          "000_OnlyLongestCDS_Erinaceus_europaeus.HEDGEHOG.cds.all.fa"],
              ["Sorex-araneus_",                    "000_OnlyLongestPEP_Sorex_araneus.COMMON_SHREW1.pep.all.fa",                                           "000_OnlyLongestCDS_Sorex_araneus.COMMON_SHREW1.cds.all.fa"],
              ["Condylura-cristata_",               "000_OnlyLongestPEP_prot_Condylura-cristata_GCF_000260355.1_ConCri1.0_rna.gbff.txt",                   "000_OnlyLongestCDS_nucl_Condylura-cristata_GCF_000260355.1_ConCri1.0_rna.gbff.txt"], 
              
              ### Chiroptera
              ["Hipposideros-armiger_",             "000_OnlyLongestPEP_prot_Hipposideros-armiger_GCF_001890085.1_ASM189008v1_rna.gbff.txt",               "000_OnlyLongestCDS_nucl_Hipposideros-armiger_GCF_001890085.1_ASM189008v1_rna.gbff.txt"], 
              ["Rhinolophus-sinicus_",              "000_OnlyLongestPEP_prot_Rhinolophus-sinicus_GCF_001888835.1_ASM188883v1_rna.gbff.txt",                "000_OnlyLongestCDS_nucl_Rhinolophus-sinicus_GCF_001888835.1_ASM188883v1_rna.gbff.txt"], 
              ["Rousettus-aegyptiacus_",            "000_OnlyLongestPEP_prot_Rousettus-aegyptiacus_GCF_001466805.2_Raegyp2.0_rna.gbff.txt",                "000_OnlyLongestCDS_nucl_Rousettus-aegyptiacus_GCF_001466805.2_Raegyp2.0_rna.gbff.txt"], 
              ["Pteropus-alecto_",                  "000_OnlyLongestPEP_prot_Pteropus-alecto_GCF_000325575.1_ASM32557v1_rna.gbff.txt",                     "000_OnlyLongestCDS_nucl_Pteropus-alecto_GCF_000325575.1_ASM32557v1_rna.gbff.txt"], 
              ["Pteropus-vampyrus_",                "000_OnlyLongestPEP_Pteropus_vampyrus.pteVam1.pep.all.fa",                                             "000_OnlyLongestCDS_Pteropus_vampyrus.pteVam1.cds.all.fa"],
              ["Desmodus-rotundus_",                "000_OnlyLongestPEP_prot_Desmodus-rotundus_GCF_002940915.1_ASM294091v2_rna.gbff.txt",                  "000_OnlyLongestCDS_nucl_Desmodus-rotundus_GCF_002940915.1_ASM294091v2_rna.gbff.txt"], 
              ["Phyllostomus-discolor_",            "000_OnlyLongestPEP_prot_Phyllostomus-discolor_GCF_004126475.1_mPhyDis1_v1.p_rna.gbff.txt",            "000_OnlyLongestCDS_nucl_Phyllostomus-discolor_GCF_004126475.1_mPhyDis1_v1.p_rna.gbff.txt"], 
              ["Miniopterus-natalensis_",           "000_OnlyLongestPEP_prot_Miniopterus-natalensis_GCF_001595765.1_Mnat.v1_rna.gbff.txt",                 "000_OnlyLongestCDS_nucl_Miniopterus-natalensis_GCF_001595765.1_Mnat.v1_rna.gbff.txt"], 
              ["Eptesicus-fuscus_",                 "000_OnlyLongestPEP_prot_Eptesicus-fuscus_GCF_000308155.1_EptFus1.0_rna.gbff.txt",                     "000_OnlyLongestCDS_nucl_Eptesicus-fuscus_GCF_000308155.1_EptFus1.0_rna.gbff.txt"], 
              ["Myotis-lucifugus_",                 "000_OnlyLongestPEP_Myotis_lucifugus.Myoluc2.0.pep.all.fa",                                            "000_OnlyLongestCDS_Myotis_lucifugus.Myoluc2.0.cds.all.fa"],
              
              ## Perissodactyla
              ["Ceratotherium-simum_",              "000_OnlyLongestPEP_Ceratotherium-simum_rna.gbk.txt",                                                  "000_OnlyLongestCDS_Ceratotherium-simum_rna.gbk.txt"],
              ["Equus-asinus-asinus_",              "000_OnlyLongestPEP_Equus_asinus_asinus.ASM303372v1.pep.all.fa",                                       "000_OnlyLongestCDS_Equus_asinus_asinus.ASM303372v1.cds.all.fa"],
              ["Equus-caballus_",                   "000_OnlyLongestPEP_Equus_caballus.EquCab2.pep.all.fa",                                                "000_OnlyLongestCDS_Equus_caballus.EquCab2.cds.all.fa"],
              ["Equus-przewalskii_",                "000_OnlyLongestPEP_prot_Equus-przewalskii_GCF_000696695.1_Burgud_rna.gbff.txt",                       "000_OnlyLongestCDS_nucl_Equus-przewalskii_GCF_000696695.1_Burgud_rna.gbff.txt"], 
              ["Manis-javanica_",                   "000_OnlyLongestPEP_Manis-javanica_rna.gbk.txt",                                                       "000_OnlyLongestCDS_Manis-javanica_rna.gbk.txt"],

              ### Carnivora
              ["Suricata-suricatta_",               "000_OnlyLongestPEP_prot_Suricata-suricatta_GCF_006229205.1_meerkat_22Aug2017_6uvM2_HiC_rna.gbff.txt", "000_OnlyLongestCDS_nucl_Suricata-suricatta_GCF_006229205.1_meerkat_22Aug2017_6uvM2_HiC_rna.gbff.txt"], 
              ["Panthera-pardus_",                  "000_OnlyLongestPEP_Panthera_pardus.PanPar1.0.pep.all.fa",                                             "000_OnlyLongestCDS_Panthera_pardus.PanPar1.0.cds.all.fa"],
              ["Panthera-pardus-1_",                "000_OnlyLongestPEP_prot_Panthera-pardus-1_GCF_001857705.1_PanPar1.0_rna.gbff.txt",                    "000_OnlyLongestCDS_nucl_Panthera-pardus-1_GCF_001857705.1_PanPar1.0_rna.gbff.txt"], 
              ["Panthera-tigris-altaica_",          "000_OnlyLongestPEP_Panthera_tigris_altaica.PanTig1.0.pep.all.fa",                                     "000_OnlyLongestCDS_Panthera_tigris_altaica.PanTig1.0.cds.all.fa"],
              ["Panthera-tigris-altaica-1_",        "000_OnlyLongestPEP_prot_Panthera-tigris-altaica-1_GCF_000464555.1_PanTig1.0_rna.gbff.txt",            "000_OnlyLongestCDS_nucl_Panthera-tigris-altaica-1_GCF_000464555.1_PanTig1.0_rna.gbff.txt"], 
              ["Lynx-canadensis_",                  "000_OnlyLongestPEP_prot_Lynx-canadensis_GCF_007474595.1_mLynCan4_v1.p_rna.gbff.txt",                  "000_OnlyLongestCDS_nucl_Lynx-canadensis_GCF_007474595.1_mLynCan4_v1.p_rna.gbff.txt"], 
              ["Acinonyx-jubatus_",                 "000_OnlyLongestPEP_prot_Acinonyx-jubatus_GCF_003709585.1_Aci_jub_2_rna.gbff.txt",                     "000_OnlyLongestCDS_nucl_Acinonyx-jubatus_GCF_003709585.1_Aci_jub_2_rna.gbff.txt"], 
              ["Puma-concolor_",                    "000_OnlyLongestPEP_prot_Puma-concolor_GCF_003327715.1_PumCon1.0_rna.gbff.txt",                        "000_OnlyLongestCDS_nucl_Puma-concolor_GCF_003327715.1_PumCon1.0_rna.gbff.txt"], 
              ["Felis-catus_",                      "000_OnlyLongestPEP_Felis_catus-E98.Felis_catus_9.0.pep.all.fa",                                       "000_OnlyLongestCDS_Felis_catus-E98.Felis_catus_9.0.cds.all.fa"],                
              ["Felis-catus-R_",                    "000_OnlyLongestPEP_Felis-catus-R_GCF_000181335.3_Felis_catus_9.0_rna.gbff.txt",                       "000_OnlyLongestCDS_Felis-catus-R_GCF_000181335.3_Felis_catus_9.0_rna.gbff.txt"],                  
              ["Vulpes-vulpes_",                    "000_OnlyLongestPEP_Vulpes_vulpes.VulVul2.2.pep.all.fa",                                               "000_OnlyLongestCDS_Vulpes_vulpes.VulVul2.2.cds.all.fa"],
              ["Canis-lupus_",                      "000_OnlyLongestPEP_prot_Canis-lupus_GCF_003254725.1_ASM325472v1_rna.gbff.txt",                        "000_OnlyLongestCDS_nucl_Canis-lupus_GCF_003254725.1_ASM325472v1_rna.gbff.txt"], 
              ["Canis-familiaris_",                 "000_OnlyLongestPEP_Canis_familiaris.CanFam3.1.pep.all.fa",                                            "000_OnlyLongestCDS_Canis_familiaris.CanFam3.1.cds.all.fa"],
              ["Enhydra-lutris_",                   "000_OnlyLongestPEP_prot_Enhydra-lutris_GCF_002288905.1_ASM228890v2_rna.gbff.txt",                     "000_OnlyLongestCDS_nucl_Enhydra-lutris_GCF_002288905.1_ASM228890v2_rna.gbff.txt"], 
              ["Neovison-vison_",                   "000_OnlyLongestPEP_Neovison_vison.NNQGG.v01.pep.all.fa",                                              "000_OnlyLongestCDS_Neovison_vison.NNQGG.v01.cds.all.fa"],
              ["Mustela-putorius-furo_",            "000_OnlyLongestPEP_Mustela_putorius_furo.MusPutFur1.0.pep.all.fa",                                    "000_OnlyLongestCDS_Mustela_putorius_furo.MusPutFur1.0.cds.all.fa"],
              ["Ailuropoda-melanoleuca_",           "000_OnlyLongestPEP_Ailuropoda_melanoleuca.ailMel1.pep.all.fa",                                        "000_OnlyLongestCDS_Ailuropoda_melanoleuca.ailMel1.cds.all.fa"],
              ["Ursus-americanus_",                 "000_OnlyLongestPEP_Ursus_americanus.ASM334442v1.pep.all.fa",                                          "000_OnlyLongestCDS_Ursus_americanus.ASM334442v1.cds.all.fa"],
              ["Ursus-arctos_",                     "000_OnlyLongestPEP_prot_Ursus-arctos_GCF_003584765.1_ASM358476v1_rna.gbff.txt",                       "000_OnlyLongestCDS_nucl_Ursus-arctos_GCF_003584765.1_ASM358476v1_rna.gbff.txt"], 
              ["Ursus-maritimus_",                  "000_OnlyLongestPEP_Ursus-maritimus_rna.gbk.txt",                                                      "000_OnlyLongestCDS_Ursus-maritimus_rna.gbk.txt"],
              ["Neomonachus-schauinslandi-1_",      "000_OnlyLongestPEP_prot_Neomonachus-schauinslandi-1_GCF_002201575.1_ASM220157v1_rna.gbff.txt",        "000_OnlyLongestCDS_nucl_Neomonachus-schauinslandi-1_GCF_002201575.1_ASM220157v1_rna.gbff.txt"], 
              ["Leptonychotes-weddellii_",          "000_OnlyLongestPEP_Leptonychotes-weddellii-R_GCF_000349705.1_LepWed1.0_rna.gbff.txt",                 "000_OnlyLongestCDS_Leptonychotes-weddellii-R_GCF_000349705.1_LepWed1.0_rna.gbff.txt"],      

              ["Odobenus-rosmarus-1_",              "000_OnlyLongestPEP_prot_Odobenus-rosmarus-1_GCF_000321225.1_Oros_1.0_rna.gbff.txt",                   "000_OnlyLongestCDS_nucl_Odobenus-rosmarus-1_GCF_000321225.1_Oros_1.0_rna.gbff.txt"], 
              ["Zalophus-californianus_",           "000_OnlyLongestPEP_prot_Zalophus-californianus_GCF_900631625.1_zalCal2.2_rna.gbff.txt",               "000_OnlyLongestCDS_nucl_Zalophus-californianus_GCF_900631625.1_zalCal2.2_rna.gbff.txt"], 
              ["Eumetopias-jubatus_",               "000_OnlyLongestPEP_prot_Eumetopias-jubatus_GCF_004028035.1_ASM402803v1_rna.gbff.txt",                 "000_OnlyLongestCDS_nucl_Eumetopias-jubatus_GCF_004028035.1_ASM402803v1_rna.gbff.txt"], 
              ["Callorhinus-ursinus_",              "000_OnlyLongestPEP_prot_Callorhinus-ursinus_GCF_003265705.1_ASM326570v1_rna.gbff.txt",                "000_OnlyLongestCDS_nucl_Callorhinus-ursinus_GCF_003265705.1_ASM326570v1_rna.gbff.txt"], 

              #### Cetartiodactyla
              ["Vicugna-pacos_",                    "000_OnlyLongestPEP_Vicugna_pacos.vicPac1.pep.all.fa",                                                 "000_OnlyLongestCDS_Vicugna_pacos.vicPac1.cds.all.fa"],
              ["Camelus-dromedarius_",              "000_OnlyLongestPEP_prot_Camelus-dromedarius_GCF_000803125.2_CamDro3_rna.gbff.txt",                    "000_OnlyLongestCDS_nucl_Camelus-dromedarius_GCF_000803125.2_CamDro3_rna.gbff.txt"], 
              ["Camelus-bactrianus_",               "000_OnlyLongestPEP_Camelus-bactrianus_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Camelus-bactrianus_rna.gbk.txt"],
              ["Sus-scrofa_",                       "000_OnlyLongestPEP_Sus_scrofa.Sscrofa11.1.pep.all.fa",                                                "000_OnlyLongestCDS_Sus_scrofa.Sscrofa11.1.cds.all.fa"],
              ["Sus-scrofa-1_",                     "000_OnlyLongestPEP_Sus-scrofa_pig.1.rna.gbff.txt",                                                    "000_OnlyLongestCDS_Sus-scrofa_pig.1.rna.gbff.txt"],
              
              #### Ruminantia
              ["Odocoileus-virginianus_",           "000_OnlyLongestPEP_Odocoileus-virginianus_rna.gbk.txt",                                               "000_OnlyLongestCDS_Odocoileus-virginianus_rna.gbk.txt"],
              ["Pantholops-hodgsonii_",             "000_OnlyLongestPEP_prot_Pantholops-hodgsonii_GCF_000400835.1_PHO1.0_rna.gbff.txt",                    "000_OnlyLongestCDS_nucl_Pantholops-hodgsonii_GCF_000400835.1_PHO1.0_rna.gbff.txt"], 
              ["Ovis-aries_",                       "000_OnlyLongestPEP_Ovis_aries.Oar_v3.1.pep.all.fa",                                                   "000_OnlyLongestCDS_Ovis_aries.Oar_v3.1.cds.all.fa"],
              ["Capra-hircus_",                     "000_OnlyLongestPEP_Capra_hircus.ARS1.pep.all.fa",                                                     "000_OnlyLongestCDS_Capra_hircus.ARS1.cds.all.fa"],
              ["Bubalus-bubalis_",                  "000_OnlyLongestPEP_Bubalus-bubalis_rna.gbk.txt",                                                      "000_OnlyLongestCDS_Bubalus-bubalis_rna.gbk.txt"],
              ["Bison-bison_",                      "000_OnlyLongestPEP_Bison-bison_rna.gbk.txt",                                                          "000_OnlyLongestCDS_Bison-bison_rna.gbk.txt"],
              ["Bos-mutus_",                        "000_OnlyLongestPEP_Bos_mutus.BosGru_v2.0.pep.all.fa",                                                 "000_OnlyLongestCDS_Bos_mutus.BosGru_v2.0.cds.all.fa"],
              ["Bos-taurus_",                       "000_OnlyLongestPEP_Bos_taurus.UMD3.1.pep.all.fa",                                                     "000_OnlyLongestCDS_Bos_taurus.UMD3.1.cds.all.fa"],
              ["Bos-taurus-1_",                     "000_OnlyLongestPEP_Bos-taurus_cow.1.rna.gbff.txt",                                                    "000_OnlyLongestCDS_Bos-taurus_cow.1.rna.gbff.txt"],
              
              ### Cetacea
              ["Balaenoptera-acutorostrata_",       "000_OnlyLongestPEP_Balaenoptera-acutorostrata_rna.gbk.txt",                                           "000_OnlyLongestCDS_Balaenoptera-acutorostrata_rna.gbk.txt"],
              ["Physeter-catodon_",                 "000_OnlyLongestPEP_Physeter-catodon_rna.gbk.txt",                                                     "000_OnlyLongestCDS_Physeter-catodon_rna.gbk.txt"],
              ["Lipotes-vexillifer_",               "000_OnlyLongestPEP_Lipotes-vexillifer_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Lipotes-vexillifer_rna.gbk.txt"],
              ["Neophocaena-asiaeorientalis_",      "000_OnlyLongestPEP_prot_Neophocaena-asiaeorientalis_GCF_003031525.1_V1_rna.gbff.txt",                 "000_OnlyLongestCDS_nucl_Neophocaena-asiaeorientalis_GCF_003031525.1_V1_rna.gbff.txt"], 
              ["Monodon-monoceros_",                "000_OnlyLongestPEP_prot_Monodon-monoceros_GCF_005190385.1_NGI_Narwhal_1_rna.gbff.txt",                "000_OnlyLongestCDS_nucl_Monodon-monoceros_GCF_005190385.1_NGI_Narwhal_1_rna.gbff.txt"], 
              ["Delphinapterus-leucas_",            "000_OnlyLongestPEP_prot_Delphinapterus-leucas_GCF_002288925.2_ASM228892v3_rna.gbff.txt",              "000_OnlyLongestCDS_nucl_Delphinapterus-leucas_GCF_002288925.2_ASM228892v3_rna.gbff.txt"], 
              ["Lagenorhynchus-obliquidens_",       "000_OnlyLongestPEP_prot_Lagenorhynchus-obliquidens_GCF_003676395.1_ASM367639v1_rna.gbff.txt",         "000_OnlyLongestCDS_nucl_Lagenorhynchus-obliquidens_GCF_003676395.1_ASM367639v1_rna.gbff.txt"], 
              ["Orcinus-orca_",                     "000_OnlyLongestPEP_Orcinus-orca_rna.gbk.txt",                                                         "000_OnlyLongestCDS_Orcinus-orca_rna.gbk.txt"],
              ["Globicephala-melas_",               "000_OnlyLongestPEP_prot_Globicephala-melas_GCF_006547405.1_ASM654740v1_rna.gbff.txt",                 "000_OnlyLongestCDS_nucl_Globicephala-melas_GCF_006547405.1_ASM654740v1_rna.gbff.txt"], 
              ["Tursiops-truncatus_",               "000_OnlyLongestPEP_Tursiops_truncatus.turTru1.pep.all.fa",                                            "000_OnlyLongestCDS_Tursiops_truncatus.turTru1.cds.all.fa"],

              #####Euarchontoglires
              ## Scandentia
              ["Tupaia-belangeri_",                 "000_OnlyLongestPEP_Tupaia_belangeri.TREESHREW.pep.all.fa",                                            "000_OnlyLongestCDS_Tupaia_belangeri.TREESHREW.cds.all.fa"],
              ["Tupaia-chinensis_",                 "000_OnlyLongestPEP_prot_Tupaia-chinensis_GCF_000334495.1_TupChi_1.0_rna.gbff.txt",                    "000_OnlyLongestCDS_nucl_Tupaia-chinensis_GCF_000334495.1_TupChi_1.0_rna.gbff.txt"], 

              ## Glires
              ["Oryctolagus-cuniculus_",            "000_OnlyLongestPEP_Oryctolagus_cuniculus.OryCun2.0.pep.all.fa",                                       "000_OnlyLongestCDS_Oryctolagus_cuniculus.OryCun2.0.cds.all.fa"],
              ["Ochotona-princeps_",                "000_OnlyLongestPEP_Ochotona_princeps.OchPri2.0-Ens.pep.all.fa",                                       "000_OnlyLongestCDS_Ochotona_princeps.OchPri2.0-Ens.cds.all.fa"],
              ["Marmota-flaviventris_",             "000_OnlyLongestPEP_prot_Marmota-flaviventris_GCF_003676075.1_ASM367607v1_rna.gbff.txt",               "000_OnlyLongestCDS_nucl_Marmota-flaviventris_GCF_003676075.1_ASM367607v1_rna.gbff.txt"], 
              ["Marmota-marmota_",                  "000_OnlyLongestPEP_prot_Marmota-marmota_GCF_001458135.1_marMar2.1_rna.gbff.txt",                      "000_OnlyLongestCDS_nucl_Marmota-marmota_GCF_001458135.1_marMar2.1_rna.gbff.txt"],
              ["Urocitellus-parryii_",              "000_OnlyLongestPEP_prot_Urocitellus-parryii_GCF_003426925.1_ASM342692v1_rna.gbff.txt",                "000_OnlyLongestCDS_nucl_Urocitellus-parryii_GCF_003426925.1_ASM342692v1_rna.gbff.txt"], 
              ["Ictidomys-tridecemlineatus_",       "000_OnlyLongestPEP_Ictidomys_tridecemlineatus.SpeTri2.0.pep.all.fa",                                  "000_OnlyLongestCDS_Ictidomys_tridecemlineatus.SpeTri2.0.cds.all.fa"],
              ["Fukomys-damarensis_",               "000_OnlyLongestPEP_Fukomys_damarensis.DMR_v1.0.pep.all.fa",                                           "000_OnlyLongestCDS_Fukomys_damarensis.DMR_v1.0.cds.all.fa"],
              ["Heterocephalus-glaber_",            "000_OnlyLongestPEP_Heterocephalus_glaber_female.HetGla_female_1.0.pep.all.fa",                        "000_OnlyLongestCDS_Heterocephalus_glaber.HetGla_female_1.0.cds.all.fa"],
              ["Octodon-degus_",                    "000_OnlyLongestPEP_Octodon_degus.OctDeg1.0.pep.all.fa",                                               "000_OnlyLongestCDS_Octodon_degus.OctDeg1.0.cds.all.fa"],
              ["Chinchilla-lanigera_",              "000_OnlyLongestPEP_Chinchilla_lanigera.ChiLan1.0.pep.all.fa",                                         "000_OnlyLongestCDS_Chinchilla_lanigera.ChiLan1.0.cds.all.fa"],
              ["Cavia-porcellus_",                  "000_OnlyLongestPEP_Cavia_porcellus.Cavpor3.0.pep.all.fa",                                             "000_OnlyLongestCDS_Cavia_porcellus.Cavpor3.0.cds.all.fa"],
              ["Cavia-aperea_",                     "000_OnlyLongestPEP_Cavia_aperea.CavAp1.0.pep.all.fa",                                                 "000_OnlyLongestCDS_Cavia_aperea.CavAp1.0.cds.all.fa"],
              ["Dipodomys-ordii_",                  "000_OnlyLongestPEP_Dipodomys_ordii.Dord_2.0.pep.all.fa",                                              "000_OnlyLongestCDS_Dipodomys_ordii.Dord_2.0.cds.all.fa"],
              ["Jaculus-jaculus_",                  "000_OnlyLongestPEP_Jaculus_jaculus.JacJac1.0.pep.all.fa",                                             "000_OnlyLongestCDS_Jaculus_jaculus.JacJac1.0.cds.all.fa"],
              ["Nannospalax-galili_",               "000_OnlyLongestPEP_Nannospalax_galili.S.galili_v1.0.pep.all.fa",                                      "000_OnlyLongestCDS_Nannospalax_galili.S.galili_v1.0.cds.all.fa"],
              ["Peromyscus-maniculatus_",           "000_OnlyLongestPEP_Peromyscus-maniculatus.Pman_1.0.pep.all.fa",                                       "000_OnlyLongestCDS_Peromyscus-maniculatus.Pman_1.0.cds.all.fa"],
              ["Cricetulus-griseus_",               "000_OnlyLongestPEP_Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.pep.all.fa",                             "000_OnlyLongestCDS_Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.cds.all.fa"],
              ["Rattus-norvegicus_",                "000_OnlyLongestPEP_Rattus_norvegicus.Rnor_6.0.pep.all.fa",                                            "000_OnlyLongestCDS_Rattus_norvegicus.Rnor_6.0.cds.all.fa"],
              ["Rattus-norvegicus-1_",              "000_OnlyLongestPEP_Rattus-norvegicus_rat.1.rna.gbff.txt",                                             "000_OnlyLongestCDS_Rattus-norvegicus_rat.1.rna.gbff.txt"],
              ["Mus-pahari_",                       "000_OnlyLongestPEP_Mus_pahari.PAHARI_EIJ_v1.1.pep.all.fa",                                            "000_OnlyLongestCDS_Mus_pahari.PAHARI_EIJ_v1.1.cds.all.fa"],
              ["Mus-spretus_",                      "000_OnlyLongestPEP_Mus_spretus.SPRET_EiJ_v1.pep.all.fa",                                              "000_OnlyLongestCDS_Mus_spretus.SPRET_EiJ_v1.cds.all.fa"],
              ["Mus-musculus_",                     "000_OnlyLongestPEP_Mus_musculus.GRCm38.pep.all.fa",                                                   "000_OnlyLongestCDS_Mus_musculus.GRCm38.cds.all.fa"],
              ["Mus-musculus-1_",                   "000_OnlyLongestPEP_Mus-musculus_mouse.rna.gbff.txt",                                                  "000_OnlyLongestCDS_Mus-musculus_mouse.rna.gbff.txt"],

              ### Dermoptera
              ["Galeopterus-variegatus_",           "000_OnlyLongestPEP_prot_Galeopterus-variegatus_GCF_000696425.1_3.0.2_rna.gbff.txt",                   "000_OnlyLongestCDS_nucl_Galeopterus-variegatus_GCF_000696425.1_3.0.2_rna.gbff.txt"], 

              ### Primates
              ["Otolemur-garnettii_",               "000_OnlyLongestPEP_Otolemur_garnettii.OtoGar3.pep.all.fa",                                            "000_OnlyLongestCDS_Otolemur_garnettii.OtoGar3.cds.all.fa"],
              ["Microcebus-murinus_",               "000_OnlyLongestPEP_Microcebus-murinus_rna.gbk.txt",                                                   "000_OnlyLongestCDS_Microcebus-murinus_rna.gbk.txt"],
              ["Propithecus-coquereli_",            "000_OnlyLongestPEP_Propithecus_coquereli.Pcoq_1.0.pep.all.fa",                                        "000_OnlyLongestCDS_Propithecus_coquereli.Pcoq_1.0.cds.all.fa"],
              ["Prolemur-simus_",                   "000_OnlyLongestPEP_Prolemur_simus.Prosim_1.0.pep.all.fa",                                             "000_OnlyLongestCDS_Prolemur_simus.Prosim_1.0.cds.all.fa"],
              ["Carlito-syrichta_",                 "000_OnlyLongestPEP_Carlito_syrichta.Tarsius_syrichta-2.0.1.pep.all.fa",                               "000_OnlyLongestCDS_Carlito_syrichta.Tarsius_syrichta-2.0.1.cds.all.fa"],
              ["Aotus-nancymaae_",                  "000_OnlyLongestPEP_Aotus_nancymaae.Anan_2.0.pep.all.fa",                                              "000_OnlyLongestCDS_Aotus_nancymaae.Anan_2.0.cds.all.fa"],
              ["Callithrix-jacchus_",               "000_OnlyLongestPEP_Callithrix_jacchus.ASM275486v1.pep.all.fa",                                        "000_OnlyLongestCDS_Callithrix_jacchus.ASM275486v1.cds.all.fa"],
              ["Cebus-capucinus_",                  "000_OnlyLongestPEP_Cebus_capucinus.Cebus_imitator-1.0.pep.all.fa",                                    "000_OnlyLongestCDS_Cebus_capucinus.Cebus_imitator-1.0.cds.all.fa"],
              ["Saimiri-boliviensis-boliviensis_",  "000_OnlyLongestPEP_Saimiri_boliviensis_boliviensis.SaiBol1.0.pep.all.fa",                             "000_OnlyLongestCDS_Saimiri_boliviensis_boliviensis.SaiBol1.0.cds.all.fa"],
              ["Colobus-angolensis-palliatus_",     "000_OnlyLongestPEP_Colobus_angolensis_palliatus.Cang.pa_1.0.pep.all.fa",                              "000_OnlyLongestCDS_Colobus_angolensis_palliatus.Cang.pa_1.0.cds.all.fa"],
              ["Rhinopithecus-roxellana_",          "000_OnlyLongestPEP_Rhinopithecus_roxellana.Rrox_v1.pep.all.fa",                                       "000_OnlyLongestCDS_Rhinopithecus_roxellana.Rrox_v1.cds.all.fa"],
              ["Rhinopithecus-bieti_",              "000_OnlyLongestPEP_Rhinopithecus_bieti.ASM169854v1.pep.all.fa",                                       "000_OnlyLongestCDS_Rhinopithecus_bieti.ASM169854v1.cds.all.fa"],
              ["Chlorocebus-sabaeus_",              "000_OnlyLongestPEP_Chlorocebus_sabaeus.ChlSab1.1.pep.all.fa",                                         "000_OnlyLongestCDS_Chlorocebus_sabaeus.ChlSab1.1.cds.all.fa"],
              ["Papio-anubis_",                     "000_OnlyLongestPEP_Papio_anubis.Panu_3.0.pep.all.fa",                                                 "000_OnlyLongestCDS_Papio_anubis.Panu_3.0.cds.all.fa"],
              ["Cercocebus-atys_",                  "000_OnlyLongestPEP_Cercocebus_atys.Caty_1.0.pep.all.fa",                                              "000_OnlyLongestCDS_Cercocebus_atys.Caty_1.0.cds.all.fa"],
              ["Mandrillus-leucophaeus_",           "000_OnlyLongestPEP_Mandrillus_leucophaeus.Mleu.le_1.0.pep.all.fa",                                    "000_OnlyLongestCDS_Mandrillus_leucophaeus.Mleu.le_1.0.cds.all.fa"],
              ["Macaca-nemestrina_",                "000_OnlyLongestPEP_Macaca_nemestrina.Mnem_1.0.pep.all.fa",                                            "000_OnlyLongestCDS_Macaca_nemestrina.Mnem_1.0.cds.all.fa"],
              ["Macaca-fascicularis_",              "000_OnlyLongestPEP_Macaca_fascicularis.Macaca_fascicularis_5.0.pep.all.fa",                           "000_OnlyLongestCDS_Macaca_fascicularis.Macaca_fascicularis_5.0.cds.all.fa"],
              ["Macaca-mulatta_",                   "000_OnlyLongestPEP_Macaca_mulatta.Mmul_8.0.1.pep.all.fa",                                             "000_OnlyLongestCDS_Macaca_mulatta.Mmul_8.0.1.cds.all.fa"],
              ["Nomascus-leucogenys_",              "000_OnlyLongestPEP_Nomascus_leucogenys.Nleu_3.0.pep.all.fa",                                          "000_OnlyLongestCDS_Nomascus_leucogenys.Nleu_3.0.cds.all.fa"],
              ["Pongo-abelii_",                     "000_OnlyLongestPEP_Pongo_abelii.PPYG2.pep.all.fa",                                                    "000_OnlyLongestCDS_Pongo_abelii.PPYG2.cds.all.fa"],
              ["Gorilla-gorilla_",                  "000_OnlyLongestPEP_Gorilla_gorilla.gorGor4.pep.all.fa",                                               "000_OnlyLongestCDS_Gorilla_gorilla.gorGor4.cds.all.fa"],
              ["Pan-troglodytes_",                  "000_OnlyLongestPEP_Pan_troglodytes.Pan_tro_3.0.pep.all.fa",                                           "000_OnlyLongestCDS_Pan_troglodytes.Pan_tro_3.0.cds.all.fa"],
              ["Pan-paniscus_",                     "000_OnlyLongestPEP_Pan_paniscus.panpan1.1.pep.all.fa",                                                "000_OnlyLongestCDS_Pan_paniscus.panpan1.1.cds.all.fa"],
              ["Homo-sapiens-E102_",                 "000_OnlyLongestPEP_Homo_sapiens-E102.GRCh38.pep.all.fa",                                             "000_OnlyLongestCDS_Homo_sapiens-E102.GRCh38.cds.all.fa"],
              #["Homo-sapiens-E99_",                 "000_OnlyLongestPEP_Homo_sapiens-E99.GRCh38.pep.all.fa",                                               "000_OnlyLongestCDS_Homo_sapiens-E99.GRCh38.cds.all.fa"],
              #["Homo-sapiens_",                     "000_OnlyLongestPEP_Homo_sapiens.GRCh38.pep.all.fa",                                                   "000_OnlyLongestCDS_Homo_sapiens.GRCh38.cds.all.fa"],
              #["Homo-sapiens-1_",                   "000_OnlyLongestPEP_Homo-sapiens_human.rna.gbff.txt",                                                  "000_OnlyLongestCDS_Homo-sapiens_human.rna.gbff.txt"],
              ["Homo-sapiens-RS2_",                 "000_OnlyLongestPEP_Homo-sapiens-RS2_GCF_000001405.40_GRCh38.p14_rna.gbff.txt",                        "000_OnlyLongestCDS_Homo-sapiens-RS2_GCF_000001405.40_GRCh38.p14_rna.gbff.txt"],

]


geneticCode = {
           #"((U|T|C|Y|A)(U|T)G)"               : "B",
           #"((C(U|T).)|((U|T)(U|T)(A|G|R)))"   : "L",
            "CTA"                             : "L",
            "CTT"                             : "L",
            "CTG"                             : "L",
            "CTC"                             : "L",
            "TTA"                             : "L",
            "TTG"                             : "L",
            "TTR"                             : "L",
           # "((CG.)|(AG(A|G|R)))"               : "R",
            "CGA"                             : "R",
            "CGT"                             : "R",
            "CGG"                             : "R",
            "CGC"                             : "R",
            "AGA"                             : "R",
            "AGG"                             : "R",
            "AGR"                             : "R",
           # "(((U|T)C.)|(AG(U|T|C|Y)))"         : "S",
            "TCA"                             : "S",
            "TCT"                             : "S",
            "TCG"                             : "S",
            "TCC"                             : "S",
            "AGT"                             : "S",
            "AGC"                             : "S",
            "AGY"                             : "S",
           # "(GC.)"                             : "A",
            "GCA"                             : "A",
            "GCT"                             : "A",
            "GCG"                             : "A",
            "GCC"                             : "A",
           # "(GG.)"                             : "G",
            "GGT"                             : "G",
            "GGC"                             : "G",
            "GGA"                             : "G",
            "GGG"                             : "G",
           # "(CC.)"                             : "P",
            "CCT"                             : "P",
            "CCC"                             : "P",
            "CCA"                             : "P",
            "CCG"                             : "P",
           # "(AC.)"                             : "T",
            "ACT"                             : "T",
            "ACC"                             : "T",
            "ACA"                             : "T",
            "ACG"                             : "T",
           # "(G(U|T).)"                         : "V",
            "GTT"                             : "V",
            "GTC"                             : "V",
            "GTA"                             : "V",
            "GTG"                             : "V",
           # "(A(U|T)(U|T|C|Y|A))"               : "I",
            "ATT"                             : "I",
            "ATC"                             : "I",
            "ATA"                             : "I",
           # "(((U|T)A(A|G|R))|((T|U)GA))"    : "_",
            "TAA"                             : "X",  #* Ter
            "TAG"                             : "X",  #* Ter
            "TAR"                             : "X",  #* Ter
            "TGA"                             : "X",  #* Ter
           # "((U|T)G(U|T|C|Y))"                 : "C",
            "TGT"                             : "C",  #Cys  
            "TGC"                             : "C",  #Cys  
            "TGY"                             : "C",  #Cys  
           # "(GA(U|T|C|Y))"                     : "D",
            "GAT"                             : "D",  #Asp
            "GAC"                             : "D",  #Asp
            "GAY"                             : "D",  #Asp
           # "(GA(A|G|R))"                       : "E",
            "GAA"                             : "E",  #Glu
            "GAG"                             : "E",  #Glu
            "GAR"                             : "E",  #Glu
           # "((U|T)(U|T)(U|T|C|Y))"             : "F",
            "TTT"                             : "F",  #Phe
            "TTC"                             : "F",  #Phe
            "TTY"                             : "F",  #Phe
           # "(CA(U|T|C|Y))"                     : "H",
            "CAT"                             : "H",  #His
            "CAC"                             : "H",  #His
            "CAY"                             : "H",  #His
           # "(AA(A|G|R))"                       : "K",
            "AAA"                             : "K",  #Lys
            "AAG"                             : "K",  #Lys
            "AAR"                             : "K",  #Lys
           # "(AA(U|T|C|Y))"                     : "N",
            "AAT"                             : "N",  #Asn
            "AAC"                             : "N",  #Asn
            "AAY"                             : "N",  #Asn
           # "(CA(A|G|R))"                       : "Q",
            "CAA"                             : "Q",  #Gln
            "CAG"                             : "Q",  #Gln
            "CAR"                             : "Q",  #Gln
           # "((U|T)A(U|T|C|Y))"                 : "Y",
            "TAT"                             : "Y",  #Tyr
            "TAC"                             : "Y",  #Tyr
            "TAY"                             : "Y",  #Tyr
           # "(A(U|T)G)"                         : "M",
            "ATG"                             : "M",  #Met
           # "((U|T)GG)"                         : "W",
            "TGG"                             : "W",  #Trp  
           # "..."                               : "X",
           # "(NNN)"                             : "X",
            "NNN"                             : "X",  
           # "(N(.|N).)"                         : "X",
            "N.."                             : "X",  
            "NN."                             : "X",  
           # "(.(.|N)N)"                         : "X",
            ".NN"                             : "X",  
            "..N"                             : "X",  
           # "(.N.)"                             : "X",
            ".N."                             : "X",  
            "---"                             : "-"}
### Codes follow pal2nal.

geneticCodeHTML = {
            'TCA' : '<span class=ser>CODON</span>',    #'S',    # ser
            'TCC' : '<span class=ser>CODON</span>',    #'S',    # ser
            'TCG' : '<span class=ser>CODON</span>',    #'S',    # ser
            'TCT' : '<span class=ser>CODON</span>',    #'S',    # ser
            'TTC' : '<span class=phe>CODON</span>',    # F
            'TTT' : '<span class=phe>CODON</span>',    #'F',    # phe
            'TTA' : '<span class=leu>CODON</span>',    #'L',    # leu
            'TTG' : '<span class=leu>CODON</span>',    #'L',    # leu
            'TAC' : '<span class=tyr>CODON</span>',    #'Y',    # tyr
            'TAT' : '<span class=tyr>CODON</span>',    #'Y',    # tyr
            'TAA' : '<span class=stp>CODON</span>',    #'*',    # stp
            'TAG' : '<span class=stp>CODON</span>',    #'*',    # stp
            'TGC' : '<span class=cys>CODON</span>',    #'C',    # cys
            'TGT' : '<span class=cys>CODON</span>',    #'C',    # cys
            'TGA' : '<span class=stp>CODON</span>',    #'*',    # stp
            'TGG' : '<span class=trp>CODON</span>',    #'W',    # trp
            'CTA' : '<span class=leu>CODON</span>',    #'L',    # leu
            'CTC' : '<span class=leu>CODON</span>',    #'L',    # leu
            'CTG' : '<span class=leu>CODON</span>',    #'L',    # leu
            'CTT' : '<span class=leu>CODON</span>',    #'L',    # leu
            'CCA' : '<span class=pro>CODON</span>',    #'P',    # pro
            'CCC' : '<span class=pro>CODON</span>',    #'P',    # pro
            'CCG' : '<span class=pro>CODON</span>',    #'P',    # pro
            'CCT' : '<span class=pro>CODON</span>',    #'P',    # pro
            'CAC' : '<span class=his>CODON</span>',    #'H',    # his
            'CAT' : '<span class=his>CODON</span>',    #'H',    # his
            'CAA' : '<span class=gln>CODON</span>',    #'Q',    # gln
            'CAG' : '<span class=gln>CODON</span>',    #'Q',    # gln
            'CGA' : '<span class=arg>CODON</span>',    #'R',    # arg
            'CGC' : '<span class=arg>CODON</span>',    #'R',    # arg
            'CGG' : '<span class=arg>CODON</span>',    #'R',    # arg
            'CGT' : '<span class=arg>CODON</span>',    #'R',    # arg
            'ATA' : '<span class=ile>CODON</span>',    #'I',    # ile
            'ATC' : '<span class=ile>CODON</span>',    #'I',    # ile
            'ATT' : '<span class=ile>CODON</span>',    #'I',    # ile
            'ATG' : '<span class=met>CODON</span>',    #'M',    # met
            'ACA' : '<span class=thr>CODON</span>',    #'T',    # thr
            'ACC' : '<span class=thr>CODON</span>',    #'T',    # thr
            'ACG' : '<span class=thr>CODON</span>',    #'T',    # thr
            'ACT' : '<span class=thr>CODON</span>',    #'T',    # thr
            'AAC' : '<span class=asn>CODON</span>',    #'N',    # asn
            'AAT' : '<span class=asn>CODON</span>',    #'N',    # asn
            'AAA' : '<span class=lys>CODON</span>',    #'K',    # lys
            'AAG' : '<span class=lys>CODON</span>',    #'K',    # lys
            'AGC' : '<span class=ser>CODON</span>',    #'S',    # ser
            'AGT' : '<span class=ser>CODON</span>',    #'S',    # ser
            'AGA' : '<span class=arg>CODON</span>',    #'R',    # arg
            'AGG' : '<span class=arg>CODON</span>',    #'R',    # arg
            'GTA' : '<span class=val>CODON</span>',    #'V',    # val
            'GTC' : '<span class=val>CODON</span>',    #'V',    # val
            'GTG' : '<span class=val>CODON</span>',    #'V',    # val
            'GTT' : '<span class=val>CODON</span>',    #'V',    # val
            'GCA' : '<span class=ala>CODON</span>',    #'A',    # ala
            'GCC' : '<span class=ala>CODON</span>',    #'A',    # ala
            'GCG' : '<span class=ala>CODON</span>',    #'A',    # ala
            'GCT' : '<span class=ala>CODON</span>',    #'A',    # ala
            'GAC' : '<span class=asp>CODON</span>',    #'D',    # asp
            'GAT' : '<span class=asp>CODON</span>',    #'D',    # asp
            'GAA' : '<span class=glu>CODON</span>',    #'E',    # glu
            'GAG' : '<span class=glu>CODON</span>',    #'E',    # glu
            'GGA' : '<span class=gly>CODON</span>',    #'G',    # gly
            'GGC' : '<span class=gly>CODON</span>',    #'G',    # gly
            'GGG' : '<span class=gly>CODON</span>',    #'G',    # gly
            'GGT' : '<span class=gly>CODON</span>',    #'G',    # gly
}

alignHTMLtopTMP = '''
<!DOCTYPE html>
<html>
<head>
    <meta http-equiv="Content-Type" content="text/html">
        <title>TITLE</title>
            <style type="text/css">
                STYLELINES
            </style>
    </head>
<body>
<pre>'''

styles_html_codons = '''
                .ala { background-color: Magenta; }
                .arg { background-color: #ff0080; }
                .asn { background-color: #0080ff; }
                .asp { background-color: #00ffff; }
                .cys { background-color: #80ff00; }
                .gln { background-color: #ffff00; }
                .glu { background-color: #ffffdf; }
                .gly { background-color: #7f7f7f; }
                .his { background-color: #ff80ff; }
                .ile { background-color: #ff80c0; }
                .leu { background-color: #80bfff; }
                .lys { background-color: #80ffff; }
                .met { background-color: #80ff80; }
                .phe { background-color: #ffff80; }
                .pro { background-color: #ff8080; }
                .ser { background-color: #bfbfbf; }
                .thr { background-color: #ffbfff; }
                .tyr { background-color: #ffbfdf; }
                .trp { background-color: #bfdfff; }
                .val { background-color: #bfffff; }
                .stp { background-color: #000000; color: white}
'''           


alignHTMLbottom = '''
</pre>
</body>
</html>'''


resHTMLlines = '''
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html lang="ja">

    <head>
        <meta http-equiv="content-type" content="text/html;charset=shift_jis">
        <!-- <meta name="generator" content="Adobe GoLive"> -->
        <title>TITLE</title>
        <link href="main.css" rel="stylesheet" type="text/css" media="all">
    </head>


 <body bgcolor="#eeeeee" leftmargin="20" marginheight="20" marginwidth="20" topmargin="20">
 <table align="center" border="0" cellspacing="0" cellpadding="5" bgcolor="white">

   <tr>
      <td width="600">
      <div align="center">
          <table width="100%" border="0" cellspacing="2" cellpadding="0" bgcolor="#000088" height="50">
      <tr>
          <td align="center" valign="middle">
          <div align="center">
              <font size="5" color=#FFFFFF face="Verdana, Arial, Helvetica, sans-serif"><b>TITLE</b></font>
          </div>
          </td>
      </tr>
      </table>
       </div>
   </td>
   </tr>

   <tr>
   <td>
     <table border="0" cellpadding="0" cellspacing="0" bordercolor="#FFFFFF">

       <tr><td align="left" valign="top">RESULTMESSAGE</td>

       </tr>
       <tr>
           <td align="left" valign="top">
             DOWNLOADMESSAGE
           </td>
       </tr>
  
       <tr>
         <td align="center" valign="top"><Hr></td>
       </tr>
  
       <!--
       <tr>
         <td align="left" valign="top"><img src="115_1stRearranged_geneTree.png"></td>
       </tr>
       -->

       <tr><td align="left" valign="top" name="REARRANGED1">REARRANGED_TREE</td></tr>
       
       <tr><td align="center" valign="top" name="REARRANGED2">Rearranged gene tree</td></tr>

       <tr><td align="center" valign="top" name="REARRANGED3"><Hr></td></tr>
  
       <tr><td align="center" valign="top" name="REARRANGED4">&#8593;</td></tr>

       <tr><td align="center" valign="top" name="REARRANGED5"><Hr></td></tr>
  
       <tr>
         <td align="left" valign="top"><img src="115_1stGeneTree.png"></td>
       </tr>
       <tr>
         <td align="center" valign="top">NJ tree</td>
       </tr>

       <tr>
         <td>
           <hr>
           <table width="100%" border="0" cellspacing="2" cellpadding="0">
             <tr>
               <td width="20%">&nbsp;</td>
               <td align="center" width="40%"><a href="http://www.fish-evol.org/index_eng.html" target="_blank"><em>jinoueATg.ecc.u-tokyo.ac.jp</em></a></td>
               <td width="20%"><div align="right"> <a href="link"></a></div></td>
             </tr>
           </table>
         </td>
       </tr> 
     </table>
   </td>
 </tr>
 </table>
    
    <br>
    <br>
</body>

</html>
'''


resHTMLlines_error = '''
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html>
    <head>
        <meta http-equiv="content-type" content="text/html;charset=shift_jis">
        <title>TITLE</title>
        <link href="main.css" rel="stylesheet" type="text/css" media="all">
    </head>


<body bgcolor="#eeeeee" leftmargin="20" marginheight="20" marginwidth="20" topmargin="20">
<table align="center" border="0" cellspacing="5" cellpadding="5" bgcolor="white">

  <!-- title -->
  <tr><td width="600"><table width="100%" border="0" cellspacing="2" cellpadding="0" bgcolor="#000088" height="50">
     <tr>
       <td align="center" valign="middle">
         <font size="5" color=#FFFFFF face="Verdana, Arial, Helvetica, sans-serif"><b>TITLE</b></font>
       </td>
     </tr>
  </table></td></tr>

  <!-- res table -->
  <tr><td width="600"><table border="1" cellpadding="0" cellspacing="0">

    <tr>
      <td align="left"  valign="top" width="15%">First query:</td>
      <td align="center"  valign="top" width="70%">FIRSTQUERY</td>
      <td align="right" valign="top" width="15%"><!-- NEXTPAGE --></td>
    </tr>

    <tr>
      <td align="left" valign="top">Result:</td>
      <td align="center" valign="top">RESULTMESSAGE</td>
      <td align="right" valign="top"><!-- PREVIOUSPAGE --></td>
    </tr>

    </table></td></tr>

    <tr><td><table width="100%" border="0" cellspacing="2" cellpadding="0">
      <tr>
        <td colspan="3"><hr></td>
      </tr>
      <tr>
        <td width="20%">&nbsp;</td>
        <td align="center" width="40%"><a href="http://www.geocities.jp/ancientfishtree/index_eng.html" target="_blank"><em>jun.inoueAToist.jp</em></a></td>
        <td width="20%"><div align="right"> <a href="link"></a></div></td>
      </tr>
    </table></td></tr> 

</table>

</body>

</html>
'''


### File uploading Start
#def FileExtensionGetAllowUpload(allow_ext_seq):
#    allow_ext_seq = ["txt", "fas", "fa"]
#    for v in allow_ext:
#        if (v == ext):
#            return True
#    return False


def check_feasibility_of_completion():
    if int(blastTopHits) == 3:
        if len(taxonSamplingList) > 120:
            print("Your analysis was stopped.<br>The number of species should be less than 120,<br>if you set the number of BLAST hits to report per genome as 3.<br>")
            print("The number of species: ", len(taxonSamplingList))
            exit()
    if int(blastTopHits) == 5:
        if len(taxonSamplingList) > 70:
            print("Error in your seqeucne file.<br>The number of species should be less than 70,<br>if you set the number of BLAST hits to report per genome as 5.<br>")
            print("The number of species: ", len(taxonSamplingList))
            exit()
    if int(blastTopHits) == 10:
        if len(taxonSamplingList) > 50:
            print("Error in your sequence file.<br>The number of species should be less than 40,<br>if you set the number of BLAST hits to report per genome as 10.<br>")
            print("The number of species: ", len(taxonSamplingList))
            exit()
    if int(blastTopHits) == 20:
        if len(taxonSamplingList) > 40:
            print("Error in your sequence file.<br>The number of species should be less than 40,<br>if you set the number of BLAST hits to report per genome as 20.<br>")
            print("The number of species: ", len(taxonSamplingList))
            exit()

def delete_dirs():
    now = datetime.date.today()
    for dir in os.listdir(dirAddress + "orthoscopeWork/"):
        if re.search("^\d", dir):
            mtime = datetime.date.fromtimestamp(int(os.path.getmtime(dirAddress + "orthoscopeWork/" + dir)))
            base, ext = os.path.splitext(dir)
            if (now - mtime).days >= 3:
                rm_command = "rm -r " + dirAddress + "orthoscopeWork/" + dir
                subprocess.call(rm_command, shell=True)


def StringGetFileExtension(string):
    string = re.split("\\\/", string)[-1]
    a      = re.split("\.", string)
    if (len(a) > 1):
        return a[-1]
    else:
        return None


def FileExtensionGetAllowUpload(allow_ext1, ext):
    for v in allow_ext1:
        if (v == ext):
            return True
    else:
        False


def uploadedFileSave(query):

    if not os.path.exists(eachDirAddress):
        os.mkdir(eachDirAddress)

    #if not (query.has_key("input_file")):
    #if not "input_file" in query:
    #    print ("Error: Select sequence file.")
    #    exit()

    #input_file = query['input_file']

    #if not (input_file.file):
    #    print("Error: Uploaded file does not contain TempraryFile object.")
    #    exit()

    #file_ext = StringGetFileExtension(input_file.filename)
    #allow_ext_seqFile = ["txt","fas"]
    #if not (FileExtensionGetAllowUpload(allow_ext_seqFile, file_ext)):
    #    print('Error in your sequence file: The extension should be ".txt" or ".fas".')
    #    exit()

    querySeqChr       = query.getvalue("querySeqChr")
    blastEvalue       = query.getvalue("blastEvalue")
    blastTopHits      = query.getvalue("blastTopHits")
    nonGapSiteRate    = float(query.getvalue("nonGapSiteRate"))
    dataset           = query.getvalue("dataset")
    if querySeqChr == "prot":
        dataset = "Amino_acid"
    RearrangementBSthreshold = query.getvalue("RearrangementBSthreshold")
    #print("RearrangementBSthreshold", RearrangementBSthreshold)
    #exit()

    taxonSamplingList = []
    input_file_tsampling = query['input_file_tsampling']
    if (input_file_tsampling.filename):
        #print("Selected input_file_tsampling")
        taxonSamplingListTMP = input_file_tsampling.file.read().decode()
        if len(taxonSamplingListTMP) > 0:
            check_tsampling_from_file(taxonSamplingListTMP)
            taxonSamplingListTMP = taxonSamplingListTMP.split("\n")
            taxonSamplingListTMP1 = []
            for a in taxonSamplingListTMP:
                if a == '' or a.startswith("#"):
                    continue
                taxonSamplingListTMP1.append(a)
            taxonSamplingList = taxonSamplingListTMP1
            
            ftTxSampling = open(eachDirAddress + "000_uploadedTxSampling.txt", "w")
            for a in taxonSamplingList:
                ftTxSampling.write(a + "\n")
            ftTxSampling.close()

        else:
            print("Error in your batch file for taxon sampling:<br>")
            print(input_file_tsampling.filename, "<br>")
            print("This file does not exist in the directory or does not contain any line.<br>")
            exit()
    else:
        #print("Were not selected input_file_tsampling")
        #taxonSamplingList = query.getlist("taxon_checkbox")
        taxonSamplingListTMP1 = query.getlist("taxon_checkbox")
        taxonSamplingListTMP2 = query.getlist("taxon_checkbox_focal")
        taxonSamplingList = taxonSamplingListTMP1 + taxonSamplingListTMP2
    #for sp in taxonSamplingList:
    #    print(sp, "|<br>")
    #exit()


    analysisGroupName = query.getvalue("analysisGroupName")
    reconciliation = query.getvalue("reconciliation")
    
    input_treeFile = ""
    '''
    if reconciliation == "treeSearchRearrange":
        if not (query.has_key("input_treeFile")):
            line_cpTree = "cp " + dirAddress + "orthoscope/examples/SpeciesTreeHypothesis.tre.txt " + eachDirAddress
            subprocess.call(line_cpTree, shell=True)
            print ("Default species tree was used.")
            #exit()
        else:
            input_treeFile = query['input_treeFile']
            file_ext = StringGetFileExtension(input_treeFile.filename)
            allow_ext_treeFile = ["tre","txt"]
            if not (FileExtensionGetAllowUpload(allow_ext_treeFile, file_ext)):
                print('Extension of tree file should be ".tre"')
                exit()
            fhtree = file(eachDirAddress + "000_uploadedTree.txt", "w")
        fhtree.write(input_treeFile.file.read())
        fhtree.close()
    '''

    #input_file_contents = input_file.file.read().decode()
    #if "−" in input_file_contents:
    #    print("Error in your uploaded fasta file.<br>")
    #    print("ORTHOSCOPE detects the N dash '-'.<br>")
    #    print("Unfortunatelly, the N dash is included in the PDF file of ORTHOSCOPE output.<br>")
    #    print("Please replace the N dashe with the hypen '-'.<br>")
    #    exit()

    #print("1111")
    text_area_seq = query.getvalue("text_area_seq")
    #print("text_area_seq", text_area_seq)
    #exit()
    text_area_seq = re.sub("[\[\]\(\):;, ]", "-", text_area_seq)
    fh = open(eachDirAddress + "000_uploaded.txt", "w")
    #fh.write(input_file_contents)
    fh.write(text_area_seq)
    fh.close()
    #print("3333")

    
    input_treeFile = ""
    if reconciliation == "treeSearchRearrange":
        #if not (query.has_key("input_treeFile")):
        if not "input_treeFile" in query:
            shutil.copy(dirAddress + "/" + version_orthoscope + "/examples/SpeciesTreeHypothesis.tre", eachDirAddress + "000_uploadedTree.txt")
        else:
            input_treeFile = query['input_treeFile']
            if not input_treeFile.filename:
                ## This is for FireFox, no tree selected.
                shutil.copy(dirAddress + "/" + version_orthoscope + "/examples/SpeciesTreeHypothesis.tre", eachDirAddress + "000_uploadedTree.txt")
            else:
                treefile_ext = StringGetFileExtension(input_treeFile.filename)
                allow_ext_treeFile = ["tre","txt"]
                if not (FileExtensionGetAllowUpload(allow_ext_treeFile, treefile_ext)):
                    print('Extension of tree file should be ".tre"')
                    exit()
                fhtree = open(eachDirAddress + "000_uploadedTree.txt", "w")
                fhtree.write(input_treeFile.file.read().decode())
                fhtree.close()
    

    #if reconciliation == "treeSearchRearrange":
    #    fhtree = file(eachDirAddress + "000_uploadedTree.txt", "w")
    #    fhtree.write(input_treeFile.file.read())
    #    fhtree.close()

    return taxonSamplingList, querySeqChr, blastEvalue, blastTopHits, nonGapSiteRate, dataset, RearrangementBSthreshold, analysisGroupName, reconciliation


def dbLinesMake():
    dbLinesFN = []
    for spNameTMP in taxonSamplingList:
        #print("spNameTMP", spNameTMP, "<br>")
        match = re.search("^([^_]+_)([^_]+)$", spNameTMP)
        spName = match.group(1)
        color = match.group(2)
        #print("spName", spName, "<br>")
        for dbLine in dbLinesTMP:
            #print("dbLine", dbLine, "<br>")
            if spName == dbLine[0]:
                #deLineNew = [dbLine[0], dbLine[1], dbLine[2], dbLine[3], color]
                deLineNew = [dbLine[0], dbLine[1], dbLine[2], color]
                dbLinesFN.append(deLineNew)
    #for dbLine in dbLinesFN:
    #    print dbLine, "<br>"
    #exit()
    return dbLinesFN
### File uploading End



### cDNA and AA file make Start
def checkUplodedFileAsFastaForamt():
    f = open(eachDirAddress + "000_uploaded.txt")
    lines_temp = list(f)
    f.close()
    
    lines = []
    for line1 in lines_temp:
        #print("line1", line1, "<br>")
        if not line1.startswith("#"):
            lines.append(line1)
    #for line2 in lines:
    #    print("line2", line2, "<br>")
    #exit()

    if not lines:
        print('Error: Cannot find your sequences in the text area.\n')
        exit()
    if not lines[0].startswith(">"):
        print('Error1 in your sequence file: In the fasta format, name lines should start with ">".\n')
        exit()

    recs_uploaded = readFasta_dict(eachDirAddress, "000_uploaded.txt")
    if int(len(recs_uploaded)) > 15:
        print("Error in your sequence file.<br>The number of queries should be less than 15.<br>")
        print("The number of queries:", len(recs_uploaded))
        exit()

    for name, seq in recs_uploaded.items():
        if re.search("\*", seq):
            print ("Error in the sequence of your fasta file.<br>")
            print ('Asterisks "*" was found.<br>')
            print ('Please remove this character.<br>')
            exit()
        if re.search("[,\/:\\\]", name):
            print ('Error in the name line of your fasta format.<br> The name line should not contain symbols such as ",()|[];:".', "<br>")
            print ("Please use only digits, letters, >, and _.")
            print ('For example, ">SpeciesName_GeneName" style is recommended.<br>')
            print ("Your name line:", name, "<br>")
            exit()
         
        if querySeqChr == "prot":
            #counter = Counter(seq)
            #print len(seq)/2
            #print (counter['A'] + counter['T'] + counter['G'] + counter['C'])
            #print "<br>"
            if len(seq)/2 < (seq.count('A') + seq.count('T') + seq.count('G') + seq.count('C')):
                print ('Error in your amino acid sequence file. Coding sequences should not be included.', "<br>")
                print ("Check the sequence of ", name, "<br>")
                exit()
        else:
            if re.search("[^ATGCNXatgcnx\- ]", seq):
                print ('Error in your coding (DNA) sequence file. <br>')
                print ('Uploaded sequences should be consist of A,T,G,C,N,or,X because you chose the "DNA" mode for your analysis.<br>')
                print ("Check the sequence of ", name, "<br>")
                exit()


def check_tsampling_from_file(tsFileContent):
    lines = tsFileContent.split("\n")
    for line in lines:
        #print("len(line)", len(line), "<br>")
        if len(line) == 0 or line.startswith("#"):
            #print("continue<br>")
            continue
        if line.startswith(">"):
            print("Error in your batch file for taxon sampling.<br>")
            print("For genome taxon sampling, probably you have uploaded a fasta file by mistake.<br>")
            print('Lines cannot start with ">".<br>')
            print("Your uploaded file contains:", line, "<br>")
            exit()
        elements = line.split("_")
        #print("elements", elements, "<br>")
        #print("len(elements)", len(elements), "<br>")
        if not len(elements) == 2:
            print("Error in your batch file for taxon sampling:<br>")
            print(line, "<br>")
            print("Your name line should be the following style:<br>")
            print("Homo-sapiens_Blue<br>")
            exit()
        speciesName_query = elements[0]
        color_query = elements[1]
        #print("speciesName_query", speciesName_query, "<br>")

        all_species_db = []        
        for dbLine in dbLinesTMP:
            all_species_db.append(dbLine[0])
        if not speciesName_query + "_" in all_species_db:
            print("Error in your batch file for taxon sampling.<br>")
            print("Your species name is not included in ORTHOSCOPE database:<br>")
            print(speciesName_query, "<br>")
            print("You can check species names by downloading our species tree file.<br>")
            exit()


#def ckeck_cDNAsequence(recsFN):
#    for name, seq in recsFN.items():
#        if re.search("[^ATGCNXatgcnx ]", seq):
#            print ("Error in your coding sequence in: ", name, "<br>")
#            print ("The sequence should be A,T,G,C,N,X. Analysis.<br>")
#            exit()


def nameChange_treeSearchRearrange(nameFN, countFN):
    if re.search(r"^>([^_]+)_", nameFN):
        name2 = re.sub("_.*$", "", nameFN)
        name3 = re.sub("^[^_]+_", "", nameFN)
        newName = name2 + "_YOURSEQ" + str(countFN) + "_" + name3
        newName = re.sub("-{2,}", "-", newName)
        newName = re.sub("\|", "", newName)
        #newName = re.sub("(.{60}).*", r"\1", newName)
        newName = re.sub("(.{" + lengthLimit_nameLine + "}).*", r"\1", newName)
        return newName
    else:
        print("Error in your sequence file: Check the name line, ", nameFN, "\n")
        exit()


def nameChange_treeSearchOnly(nameFN, countFN):
    newName = "YOURSEQ" +str(countFN) + "_" + nameFN[1:]
    newName = re.sub("[,\/ ]", "-", newName)
    newName = re.sub("-{2,}", "-", newName)
    newName = re.sub("\|", "", newName)
    newName = re.sub("-_", "_", newName)
    newName = re.sub("_-", "_", newName)
    #newName = re.sub("(.{60}).*", r"\1", newName)
    newName = re.sub("(.{" + lengthLimit_nameLine + "}).*", r"\1", newName)
    newName = ">" + newName
    return newName


def cDNAqueryFileMaker():

    recs_uploaded = readFasta_dict(eachDirAddress, "000_uploaded.txt")
    #ckeck_cDNAsequence(recs_uploaded)
    f = open(eachDirAddress + "000_cDNAQuery.txt", "w")

    count = 1
    for name, seq in recs_uploaded.items():
        name = name.rstrip("\n")
        name = name.rstrip("\r")
        if reconciliation == "treeSearchRearrange":
            newName2 = nameChange_treeSearchRearrange(name, count)
            f.write(newName2 + "\n")
        else:
            newName2 = nameChange_treeSearchOnly(name, count)
            f.write(newName2 + "\n")
        
        seq = re.sub("-", "", seq)
        f.write(seq + "\n")

        count += 1

    f.close()


def read_summary():
    f = open(eachDirAddress + "100_1stAnalysisSummary.txt", "r")
    InfileLines = list(f)
    sumDictFN  = OrderedDict()
    flag = 0
    for Line in InfileLines:
        Line = Line.rstrip("\n")
        if not Line:
            continue
        if Line[0] == "#":
            continue
        if Line[0] == ">":
            Name            = Line
            sumDictFN[Name] = ""
        else:
            sumDictFN[Name] += Line + "|"
        if Line == "</pre>":
            break

    return sumDictFN


def readPhy_dict(phyFileName):
    phyFile = open(eachDirAddress + phyFileName, "r")
    lines = list(phyFile)
    seqDictFN  = OrderedDict()
    for line in lines[1:]:
        name,seq = re.split(" +", line)
        seq = seq.rstrip("\n")
        seqDictFN[">" + name] = seq
    phyFile.close()
    return seqDictFN


def readFasta_dict(dirAddressFN, InfileNameFN):
    Infile = open(dirAddressFN + InfileNameFN, "r")
    seqDictFN  = OrderedDict()
    for Line in Infile:
        Line = Line.rstrip("\n")
        if not Line:
            continue
        elif Line[0] == "#":
            continue
        elif Line[0] == ">":
            Name = Line
            seqDictFN[Name] = ""
        else:
            Line = Line.replace("\n", "")
            Line = Line.replace("\r", "")
            #if InfileNameFN == "000_uploaded.txt":
            #    Line = re.sub("-", "", Line)
            seqDictFN[Name] += Line.upper()
    Infile.close()
    #for name, sec in seqDictFN.items():
    #    print("name", name, "<br>")
    #    print("sec", sec, "<br>")
    #exit()
    return seqDictFN


def splitDna(dna):
    codons = []
    for start in range(0, len(dna)-2, 3):
        codons.append(dna[start:start+3])
    return(codons)


def translation(dna):
    dna = dna.upper()
    protein = ""
    for codon in splitDna(dna):
        #aa = ""
        #selectedAA = []
        #for key in geneticCode:
        #    if re.match(key, codon):
        #        selectedAA.append(geneticCode.get(key))
        ## B or M
        #if len(selectedAA) > 1:
        #    aa = selectedAA[1]
        #else:
        #    aa = selectedAA[0]
        
        aa = geneticCode.get(codon, "X")

        protein = protein + aa
    return protein


def aaSeqMaker():
    if querySeqChr == "nucl":
        recfn = readFasta_dict(eachDirAddress, "000_cDNAQuery.txt")
    else:
        recfn = readFasta_dict(eachDirAddress, "000_uploaded.txt")

    fa = open(eachDirAddress + "000_aaQuery.txt","w")

    # >Human_ENSP00000259365 gene:ENSG00000136842 transcript:ENST00000259365 gene_biotype:protein_coding
    
    nums_querySeqs = len(recfn)

    if querySeqChr == "nucl":
        count = 1
        for name, seq in recfn.items():
            newName = re.sub("\|", "", name)
            fa.write(newName + "\n")
            fa.write(translation(seq) + "\n")
            count += 1
    else:
        count = 1
        for name, seq in recfn.items():
            newName = re.sub("\|", "", name)
            if reconciliation == "treeSearchRearrange":
                newName2 = nameChange_treeSearchRearrange(newName, count)
                fa.write(newName2 + "\n")
            else:
                newName2 = nameChange_treeSearchOnly(name, count)
                fa.write(newName2 + "\n")
             
            seq = re.sub("-", "", seq)
            fa.write(seq + "\n")
            count += 1

    fa.close()
    
    return nums_querySeqs


### cDNA and AA file make End

### Tree file

def get_nodeNames_ancestral():
    nodeNames_ancestral = ""
    if analysisGroupName == "Protostomia" or analysisGroupName == "Deuterostomia": 
        nodeNames_ancestral = ["Nephrozoa", "Bilateria"]
#        nodeNames_ancestral = ["Bilateria", "CnidariaHomo", "TrichoplaxHomo", "AmphimedonHomo", "Metazoa"]
    elif analysisGroupName == "Acropora":
        nodeNames_ancestral = ["Scleractinia","Anthozoa","Cnidaria"]
    elif analysisGroupName == "Chordata":
        #nodeNames_ancestral = ["Deuterostomia","Nephrozoa","Bilateria"]
        nodeNames_ancestral = ["Nephrozoa", "Bilateria"]
    elif analysisGroupName == "Vertebrata":
        #nodeNames_ancestral = ["Olfactore", "Chordata", "Deuterostomia", "Bilateria"]
        nodeNames_ancestral = ["Olfactore", "Chordata", "Deuterostomia"]
        #nodeNames_ancestral = ["Olfactore", "Chordata"]
    elif analysisGroupName == "Mammalia":
        nodeNames_ancestral = ["Amniota", "Tetrapoda", "Sarcopterygii"]
    elif analysisGroupName == "Plants":
        nodeNames_ancestral = ["Magnoliopsida","AmborellalesMesangiospermae","Mesangiospermae"]
    else:
        nodeNames_ancestral = ["Teleostomi"]
        nodeNames_ancestral = ["Teleostomi", "Gnathostomata", "Vertebrata"]
        #nodeNames_ancestral = ["Teleostomi", "Gnathostomata", "Vertebrata", "Olfactore", "Chordata"]
    return nodeNames_ancestral


def leafCollectInOrderFrom_bothNHXnewick (treeFileName):
    treeTMP = open(eachDirAddress + treeFileName, "r")
    tree = list(treeTMP)[0]
    treeTMP.close()
    leaves = tree.split(",")
    leaves = [re.sub("^\(*([^:\(\)]+).*$", r"\1", leaf) for leaf in leaves]
    leaves = [re.sub("[\t\n]", "", leaf) for leaf in leaves]
    #for leaf in leaves:
    #    print(leaf, "<br>")
    #exit()
    return leaves


def TEST_querySpIsInTree():
    recsUploadedCDNA = readFasta_dict(eachDirAddress, "000_uploaded.txt")
    speciesNamePrefixes = []
    
    for name, seq in recsUploadedCDNA.items():
        match = re.search(r"^>([^_ ]+)_", name)
        if not match:
            print('Error in your sequence file.<br>The name line should be ">Homo-sapiens_Brachyury." style.<br>')
            print('Use - "hyphen" between Homo and sapiens,<br>')
            print(' and _ "under bar" between Homo-sapiens and Brachyury.<br>') 
            print('Otherwise, use Estimating gene tree mode if this is your first trial.<br>')
            exit()

    for name, seq in recsUploadedCDNA.items():
        match = re.search(r"^>([^_ ]+)_", name)
        if match:
            speciesNamePrefixes.append(match.group(1))

    if not speciesNamePrefixes:
        print ("Error in your sequence file: The name line of sequence record should start with >Xxxx-xxxx_.")
        exit()

    spTreeTMP = open(eachDirAddress + "000_uploadedTree.txt")
    spTree = list(spTreeTMP)[0]
    spTreeTMP.close()
    if re.search("\)\)", spTree) or re.search("\),", spTree) :
        print ("Error in your tree file: In the newick format of ORTHOSCOPE, all nodes should have their node names. <br>")
        print ("Check you newick format using FigTree or TreeGraph_2. <br>")
        exit()
    num_rightOpen = spTree.count("(")
    num_leftOpen  = spTree.count(")")
    num_comma     = spTree.count(",")
    if num_rightOpen != num_leftOpen:
        print ("Error in your tree file: In the newick format of ORTHOSCOPE, numbers of right and left paretheses should be same. <br>")
        exit()
    if num_rightOpen != num_comma or num_leftOpen != num_comma :
        print ("Error in your tree file: In the newick format of ORTHOSCOPE, all nodes should be bifurcated. <br>")
        exit()

    allNodes  = nodesCollectFromNewick(eachDirAddress, treeFileName = "000_uploadedTree.txt")

    orthoGroupSpNames = ""
    #nodeNames_ancestral = ["Teleosto"]
    for nodeName_ancestral in nodeNames_ancestral:
        #print("nodeName_ancestral", nodeName_ancestral, "<br>")
        for node in allNodes:
            if nodeName_ancestral == node[2]:
                orthoGroupSpNames = node[1]
    if not orthoGroupSpNames:
        print ("Error: Orthogroup node, ", nodeNames_ancestral, " is not found in species tree.")
        exit()

    #print("orthoGroupSpNames", orthoGroupSpNames)
    #exit()
    if not speciesNamePrefixes[0] in orthoGroupSpNames:
        print ("Error in your query file. <br>Species name of your first query, ", speciesNamePrefixes[0], ", is not found in the ", nodeNames_ancestral[-1], " clade of the species tree.<br>")
        print ('Try the "Estimating gene tree" mode at first, if you are a new user.<br>')
        print ('Otherwise, if you selected the "Comparing gene and species trees" mode, name lines should be the "Homo-sapiens_ID_Description" style. "Homo-sapiens" is used for comparisons between gene and species trees.<br>')
       
        exit()

    allLeaves = allNodes[0][1]
    for speciesNamePrefix in speciesNamePrefixes:
        if not speciesNamePrefix in allLeaves:
            print ("Error in your seqeunce file. ", speciesNamePrefix, ", is not found in the species tree.<br>")
            print ("Add ", speciesNamePrefix, "to the species tree.<br>")
            print ('Otherwise, try "Estimating gene tree" mode.')
            print ('Please send an email to jinoue@g.ecc.u-tokyo.ac.jp.')
            exit()
    for dbLine in dbLines:
        spName_dbLine = dbLine[0]
        spName_dbLine = re.sub("_$","",spName_dbLine)
        if not spName_dbLine in allLeaves:
            print ("In the species tree, ", spName_dbLine, " is not found.")
            print ("Add this species to the species tree, or try Estimating gene tree mode.<br>")
            print ('Please send an email to jinoue@g.ecc.u-tokyo.ac.jp.')
            exit()
    return allNodes


def nodesCollectFromNewick(eachDirAddress, treeFileName):

    t = open(eachDirAddress + treeFileName)
    treeFN = "".join(list(t))
    treeFN = re.sub("[ \n]", "", treeFN)
    t.close()
    
    nodes = []     # 2D array for nodes
    cladeReg = "\(([^\(\)]+)\)([\d\w]+)"
    while re.search(cladeReg, treeFN):
        treeFN = treeFN.rstrip("\n")
        match = re.search(cladeReg, treeFN)               # Pick up the smallest and leftmost clade for the following analysis
        leavesString = match.group(1)
        exp          = match.group(2)
        leaves = set(leavesString.split(","))
        leaves = [re.sub(":.+", "", leaf) for leaf in leaves]  # Delete branch lengths
        leaves = [re.sub(" ", "", leaf) for leaf in leaves]  # Delete space
        treeFN = re.sub(cladeReg, r"\1", treeFN, 1)            # Delete the outer parentheses from the analyzed clade
        nodes.append([len(leaves), leaves, exp])

    sortedNodes = sorted(nodes, key=lambda x:x[0], reverse=True)

    # Add leaves as nodes
    largestNode = sortedNodes[0]
    for leaf in largestNode[1]:
        #tempLeafNode = {leaf}
        #tempLeafNode = set(leaf)
        tempLeafNode = [leaf]
        sortedNodes.append([1, tempLeafNode, leaf])

    return sortedNodes


def nhx2NewickWithNHXnodeName(tree):
  
    cladeReg = "\)([^\[]+\[&&NHX[^\]]+\])"
    count = 0
    while re.search(cladeReg, tree):
        #print("count   :", count)
        matchA = re.search(cladeReg, tree)
        expTemp  = matchA.group(1);
        #print("expTemp:", expTemp)
        #if count == 10:
        #    exit()

        matchB = re.search("(\[.*\])", expTemp)
        exp  = matchB.group(1)
 
        if re.search("B=([\d]+)", expTemp):
            matchC = re.search("B=([\d]+)", expTemp)
            bs = matchC.group(1)
        else:
            bs = "r"

        exp = re.sub("&&NHX:", "", exp)
        exp = re.sub(":",      "_", exp)
        exp = re.sub("_B=.*$", "", exp)
        
        #print("expTemp :", expTemp)
        #print("bs      :", bs)
        #print("exp     :", exp)
        #print("bs + exp:", bs + exp)

        tree = re.sub(cladeReg, ")" + bs + "_" + exp, tree, 1)
        count += 1
        #print("tree2: ", tree)
        #print()
    
    tree = re.sub("\:[\d|\.E-]+", "", tree)   #Delete the branch lengths

    ## left node name
    tree = re.sub("\[&&NHX:[^]]+\]", "",  tree)
    tree = re.sub("[\[\]]",          "", tree)
    #tree = re.sub("\]",              "",  tree)

    return tree;
### 


### Ortholog select Start
def cladesCollectFromNewick(treeFN):
    clades = []     # 2D array for clades
    cladeReg = "\(([^\(\)]+)\)(\w+)"
    while re.search(cladeReg, treeFN):
        treeFN       = treeFN.rstrip("\n")
        match        = re.search(cladeReg, treeFN)                # Pick up the smallest and leftmost clade for the following analysis
        leavesString = match.group(1)
        exp          = match.group(2)
        leaves       = set(leavesString.split(","))
        treeFN       = re.sub(cladeReg, r"\1", treeFN, 1)      # Delete the outer parentheses from the analyzed clade
        clades.append([len(leaves), leaves, exp])

    sortedClades     = sorted(clades, key=lambda x:x[0], reverse=True)

    ## Add leaf as a clade
    largestClade = sortedClades[0];
    for leaf in largestClade[1]:
        exp           = leaf
        tempLeafClade = set([leaf])
        sortedClades.append([1, tempLeafClade, exp])

    return sortedClades


def cladesCollectFromNHX(treeFN):

    clades = []     # 2D array for clades
    expReg = "\[.*?\]"
    cladeReg = "\(([^\(\)]+)\)(.*?\[.*?\])"
    while re.search(cladeReg, treeFN):
        treeFN       = treeFN.rstrip("\n")
        match        = re.search(cladeReg, treeFN)                # Pick up the smallest and leftmost clade for the following analysis
        leavesString = match.group(1)
        exp          = match.group(2)
        leavesString = re.sub("\[.*?\]",        "",    leavesString)   # Delete exp [...] of internal branches
        leavesString = re.sub(":\d+\.\d+E-\d*", "",    leavesString)   # Delete blanch lengths with E-
        leavesString = re.sub(":\d+\.\d+",      "",    leavesString)   # Delete blanch lengths of internal branches
        leaves       = set(leavesString.split(","))
        treeFN       = re.sub(cladeReg,    r"\1", treeFN, 1)      # Delete the outer parentheses from the analyzed clade
        clades.append([len(leaves), leaves, exp])

    sortedClades     = sorted(clades, key=lambda x:x[0], reverse=True)

    ## Add leaf as a clade
    largestClade = sortedClades[0];
    for leaf in largestClade[1]:
        #print ("leaf:", leaf)
        match         = re.search(r"([^_]+)_", leaf)
        exp           = "[&&NHX:S=" + match.group(1) +"]"
        tempLeafClade = set([leaf])
        sortedClades.append([1, tempLeafClade, exp])

    return sortedClades


def orthoGroupIdentifyFromNHX (treeNHX):

    #f1stTree = open(eachDirAddress + f1stTree)
    #treeNHX = list(f1stTree)[0]
    #f1stTree.close()

    f1stTree = open(eachDirAddress + "000_uploadedTree.txt")
    speciesTree = "".join(list(f1stTree))
    speciesTree = re.sub("[ \n]", "", speciesTree)
    f1stTree.close()

    topHits = topHitPicker()
    #for name, value in topHits.items():
    #    print name, "<br>"
    #    print value, "<br>"
    #exit()

    topHitName_1stQuery = list(topHits.values())[0][0]
    topHitName_1stQueryUploaded = list(topHits.keys())[0]
    
    #print("topHitName_1stQuery:", topHitName_1stQuery)
    #print("topHitName_1stQueryUploaded:", topHitName_1stQueryUploaded)
    #exit()

    allCladesSR   = cladesCollectFromNHX(treeNHX)


    keyClades   = []
    #print("cladeNamesSR", cladeNamesSR, "<br>")
    for clade in allCladesSR:
        #print "clade[1]:", clade[1], "<br>"
        #for leaf in clade[1]:
        #    print leaf
        #exit()
        for nodeName_ancestral in nodeNames_ancestral:
            #print("nodeName_ancestral:", nodeName_ancestral)
            if re.search("S=" + nodeName_ancestral + ":", clade[2]):
                #print("  FOUND")
                if [leaf for leaf in clade[1] if re.search(topHitName_1stQuery[1:], leaf)]:
                    keyClades.append(clade)
            #print("<br><br>")
    orthoGroupSR = []
    if not keyClades:
        if re.search("NoBlastHit", topHitName_1stQuery):
            stringTMP = "No blast hit was found for your 1st query within genome of this species."
        else:
            stringTMP = "No orthogroup was found for :" + topHitName_1stQueryUploaded + "<br> in the rearranged tree. "
            if analysisGroupName == "Protostomia":
                stringTMP += "Include deuterostomes.<br>"
            elif analysisGroupName == "Deuterostomia":
                stringTMP += "Include protostomes.<br>"
            elif analysisGroupName == "Acropora":
                stringTMP += "Include non-acroporid cnidarians.<br>"
            elif analysisGroupName == "Vertebrata":
                stringTMP += "Include non-vertebrate chordates such as Branchiostoma belcheri (lancelet).<br>"
            elif analysisGroupName == "Plants":
                stringTMP += "Include non-asterid Mesangiospermae plants.<br>"
            else:
                stringTMP += "Include tetrapods.<br>"
        orthoGroupSR = [0, 0, stringTMP]
        orthogroupNodeInfor = ""
    else:
        orthoGroupSR = keyClades.pop()
        #print("orthoGroupSR", orthoGroupSR)
    return orthoGroupSR


def isIndependent(checkCladeLeave, focalCladeLeave):
    for ckeckCladeLeaf in checkCladeLeave:
        if [focalCladeLeaf for focalCladeLeaf in focalCladeLeave if re.search(ckeckCladeLeaf, focalCladeLeaf)]:
            return 0;
    return 1;

def childCladesCollect(allCladesSR, orthoGroupSR):
    childCladesSR = []
    for eachClade in allCladesSR:
        if eachClade[1].issubset(orthoGroupSR[1]):
            childCladesSR.append(eachClade)
    return childCladesSR;


def sisterCladeIdentify(childCladesSR, daughterGroup1stSR):
    daughterGroup2ndSR = ""
    for eachChildClade in childCladesSR:
        if isIndependent(eachChildClade[1], daughterGroup1stSR[1]):
            daughterGroup2ndSR = eachChildClade;
            break
    return daughterGroup2ndSR


def daughterGroupsIdentify(allCladesSR, orthoGroupSR):
    childClades      = childCladesCollect(allCladesSR, orthoGroupSR)
    daughterGroup1st = childClades[1];
    daughterGroup2nd = sisterCladeIdentify(childClades, daughterGroup1st)
    return(daughterGroup1st, daughterGroup2nd)


def selectOrthogroupsFromNewick(treeNewick, cladeNamesSR, topHitName_1stQuery):
    keyClades   = []
    allCladesSR  = cladesCollectFromNewick (treeNewick)
    #print("cladeNamesSR", cladeNamesSR, "<br>")
    for clade in allCladesSR:
        #print "clade[2]:", clade[2], "<br>"
        for cladeNameSR in cladeNamesSR:
            if clade[2] == cladeNameSR:
                if [leaf for leaf in clade[1] if re.search(topHitName_1stQuery, leaf)]:
                    keyClades.append(clade)
    if not keyClades:
        print ("Error: No orthoGroup was found for ", cladeNamesSR, topHitName_1stQuery)
        exit()
    else:
        orthoGroupSR = keyClades.pop()
        return orthoGroupSR
### Ortholog select End


### makeSummary Start
def count_blastHits(infile):
    f = open(eachDirAddress + infile)
    lines = list(f)
    f.close()
    
    spNamePrefixes = []
    for spNameTMP in taxonSamplingList:
        match = re.search("^([^_]+_)([^_]+)$", spNameTMP)
        spName = match.group(1)
        spNamePrefixes.append(spName)

    rec_blastHitsNums = OrderedDict()
    for spNamePrefix in spNamePrefixes:
        hits = len([line for line in lines if re.search(spNamePrefix, line)])
        rec_blastHitsNums[spNamePrefix[:-1]]= hits
    return rec_blastHitsNums


def make_cladeSpNamePrefixes(speciesTree, nodeNames, topHitName_1stQuery):
    spNamePrefixesInClade = []
    spNamePrefix_topHitName_1stQuery = re.sub("^>([^_]+)_.*$", r"\1", topHitName_1stQuery)
    spNamePrefixes_clade2D = selectOrthogroupsFromNewick(speciesTree, nodeNames, spNamePrefix_topHitName_1stQuery)
    spNamePrefixesTMP = spNamePrefixes_clade2D[1]
    #print("spNamePrefixesTMP", spNamePrefixesTMP)
    #exit()
    for dbLine in dbLines:
        if dbLine[0][:-1] in spNamePrefixesTMP:
            spNamePrefixesInClade.append(dbLine[0])
    return spNamePrefixesInClade


def leafCountFromSpNamePrefix(orthoGroup, spNamePrefixes):
    rec_orthoNums = OrderedDict()
    for spNamePrefix in spNamePrefixes:
        #hits = len([orthoLeaf for orthoLeaf in orthoGroup[1] if re.search(spNamePrefix, orthoLeaf)])
        hits = 0
        for orthoLeaf in orthoGroup[1]:
            if re.search(spNamePrefix, orthoLeaf):
                hits += 1
        rec_orthoNums[spNamePrefix[:-1]]= hits
    return rec_orthoNums


def make_Res_orthogroupSuppor(orthogroupNodeInfor):
    if orthogroupNodeInfor.startswith("r"):
        orthogroupNodeInfor = "Orthogroup was delineated by the rearranged node.\n"
        return orthogroupNodeInfor
    elif re.search("B=([\d\.]+)\]", orthogroupNodeInfor):
        match = re.search("B=([\d\.]+)\]", orthogroupNodeInfor)
        bsValue = match.group(1)
        bsValue = re.sub("\.\d+", "", bsValue)
        return "The orthogroup monophyly is supported by " + str(bsValue) + "% bootstrap value." 
    else:
        return orthogroupNodeInfor


def error_makeSummary(resultFN):
    fs = open(eachDirAddress + "100_1stAnalysisSummary.txt", "w")

    fs.write("################ Results ################\n\n")
    fs.write(">Result\n")
    fs.write(resultFN + "\n");
    fs.write("\n")

    uploadedSEQfn = readFasta_dict(eachDirAddress, "000_uploaded.txt")
    fs.write(">Querys_used_in_the_analysis\n")
    lines_query = make_querySeqLines4summary(uploadedSEQfn)
    for line in lines_query:
        fs.write(line)
    fs.write("\n")

    fs.close()


def change_txt2html():
    summaryFN = read_summary()
    
    out = open(eachDirAddress + "100_1stAnalysisSummary.html", "w")
    alignHTMLtop = make_alignHTMLtop()
    out.write(alignHTMLtop + "\n")
    for name, contentsTMP in summaryFN.items():

        if name.startswith(">Number_of_orthogroup_member") or name.startswith(">Number_of_blastHits"):
            out.write(name + "\n")
            contents = contentsTMP.split("|")
            for content in contents:
                if not content:
                    continue
                speciesName = ""
                match = re.search("^([^ ]+)", content)
                speciesName = match.group(1)
                out.write("<span class=" + speciesName + "_>" + content + "</span>\n")
        elif name.startswith(">taxonSampling_color"):
            out.write(name + "\n")
            contents = contentsTMP.split("|")
            for content in contents:
                if not content:
                    continue
                speciesName = re.sub("_[^_]+$", "", content)
                out.write("<span class=" + speciesName + "_>" + content + "</span>\n")
        else:
            if name.startswith(">Queries_used_in_the_analysis"):
                out.write("################ Results ################\n\n")
            if name.startswith(">Mode"):
                out.write("################ Settings ################\n\n")
            out.write(name + "\n")
            contents = contentsTMP.split("|")
            for content in contents:
                out.write(content + "\n")
        out.write("\n")

    out.write(alignHTMLbottom + "\n")
    out.close()


def makeSummary():

    uploadedSEQfn = readFasta_dict(eachDirAddress, "000_uploaded.txt")
    topHits = topHitPicker()

    cDNAfn = ""
    if querySeqChr == "nucl":
        cDNAfn               = readFasta_dict(eachDirAddress, "000_cDNAQuery.txt")
    recAAfn              = readFasta_dict(eachDirAddress, "000_aaQuery.txt")
    rec044_unambSiteRate = readFasta_dict(eachDirAddress, "044_unambSiteRate.txt")

    fs = open(eachDirAddress + "100_1stAnalysisSummary.txt", "w")

    fs.write("################ Results ################\n\n")

    fs.write(">Queries_used_in_the_analysis\n")
    fs.write(make_line_for_topHits(uploadedSEQfn, topHits))
    fs.write("\n")

    speciesTree = ""
    if reconciliation == "treeSearchRearrange":

        f1stTree = open(eachDirAddress + "085_NJBS1st.txt.rearrange.0")
        treeNHX = list(f1stTree)[0]
        f1stTree.close()

        f1stTree = open(eachDirAddress + "000_uploadedTree.txt")
        speciesTree = "".join(list(f1stTree))
        speciesTree = re.sub("[ \n]", "", speciesTree)
        f1stTree.close()

        allClades   = cladesCollectFromNHX(treeNHX)
        orthoGroup = orthoGroupIdentifyFromNHX(treeNHX)
        orthogroupNodeInfor = ""
        if orthoGroup[0] == 0:
            fs.write(">Orthogroup:\n")
            orthoGroupTMP = orthoGroup[2]
            orthoGroupTMP = re.sub("<br>", "\n", orthoGroupTMP)
            fs.write(orthoGroupTMP + "\n")
        else:
            fs.write(">BootstrapValue_OrthogroupMonophyly\n")
            Res_orthogroupSuppor = make_Res_orthogroupSuppor(orthoGroup[2])
            fs.write(Res_orthogroupSuppor + "\n")
            fs.write("\n")

            fs.write(">Orthogroup\n")
            #for dbLine in dbLines:
            #    hitLeaves = [leaf for leaf in orthoGroup[1]if re.search(dbLine[0], leaf)]
            #    for leaf in hitLeaves:
            #        fs.write(leaf + "\n")
            orthoGroupOrderd = sorted(orthoGroup[1])
            for leaf in orthoGroupOrderd:
                fs.write(leaf + "\n")
            fs.write("\n")
            
            match = re.search("S=([^:]+):", orthoGroup[2])
            nodeName_orthogroup = match.group(1)
            #print nodeName_orthogroup
            #exit()
            topHits = topHitPicker()
            topHitName_1stQuery = list(topHits.values())[0][0]
            spNamePrefixes_orthoGroup  = make_cladeSpNamePrefixes(speciesTree, [nodeName_orthogroup], topHitName_1stQuery)
            rec_orthoNums = leafCountFromSpNamePrefix(orthoGroup, spNamePrefixes_orthoGroup)
            rec_orthoNums = whiteSpaceAdd(rec_orthoNums)
            fs.write(">Number_of_orthogroup_member\n")
            for name, num in rec_orthoNums.items():
                fs.write(name + ":" + str(num) + "\n")
            fs.write("\n")

            focalGroup = ""
            sisterGroup = ""
            daughterGroup1st, daughterGroup2nd = daughterGroupsIdentify(allClades, orthoGroup)
            if topHitName_1stQuery[1:] in daughterGroup1st[1]:
                focalGroup  = daughterGroup1st
                sisterGroup = daughterGroup2nd
            else:
                focalGroup  = daughterGroup2nd
                sisterGroup = daughterGroup1st

        fs.write(">Rearranged_gene_tree_newick\n")
        rearrangedTreeNewick = nhx2NewickWithNHXnodeName(treeNHX)
        fs.write(rearrangedTreeNewick)
        fs.write("\n")

        fs.write(">Rearranged_gene_tree_nhx\n")
        fs.write(treeNHX)
        fs.write("\n")

    fs.write(">Gene_tree_newick\n")
    f1stTree = open(eachDirAddress + "085_NJBS1st.txt")
    fs.write(list(f1stTree)[0])
    f1stTree.close()
    fs.write("\n")


    fs.write(">Number_of_blastHits\n")
    rec_blastHits = count_blastHits(infile = "010_blastRes.txt")
    rec_blastHits = whiteSpaceAdd(rec_blastHits)
    for name, num in rec_blastHits.items():
        fs.write(name + ":" + str(num) + "\n")
    fs.write("\n")

    rec044_unambSiteRate = whiteSpaceAdd(rec044_unambSiteRate)
    fs.write(">Aligned-site_rate\n")
    for name, siteRate in rec044_unambSiteRate.items():
        nameRate = name + ":" + str(siteRate)
        if float(siteRate) > nonGapSiteRate:
            fs.write(nameRate + "\n")
        else:
            fs.write("== Removed ==> " + nameRate + "\n")
    fs.write("\n")

    fs.write(">AnalysisTime\n")
    elapsed_time = round((time.time() - startTime),1)
    fs.write (str(elapsed_time) + " seconds")
    fs.write("\n\n")

    fs.write("\n################ Settings ################\n\n")

    fs.write(">Mode\n")
    if reconciliation == "treeSearchRearrange":
        fs.write("Comparing gene and species trees\n")
    else:
        fs.write("Estimating gene tree\n")
    fs.write("\n")

    fs.write(">AnalysisGroup\n")
    fs.write(analysisGroupName + "\n")
    fs.write("\n")

    fs.write(">Uploaded_your_sequences\n")
    for name, seq in uploadedSEQfn.items():
        name = re.sub("[\n\r]", "", name)
        fs.write(name[1:] + "\n")
        fs.write(seq + "\n")
    fs.write("\n")

    if reconciliation == "treeSearchRearrange":
        fs.write(">SpeciesTree\n")
        fs.write(speciesTree)
    fs.write("\n")

    fs.write(">taxonSampling_color\n")
    for spNameTMP in taxonSamplingList:
        #match = re.search("^([^_]+_)([^_]+)$", spNameTMP)
        #spName = match.group(1)
        #color = match.group(2)
        #spaces = " " * (40 - len(spName))
        #if color == "#660099":
        #    color = "Purple"
        #if color == "#FF00FF":
        #    color = "Magenta"
        #fs.write(spName + spaces + color +"\n")
        fs.write(spNameTMP + "\n")
    #for dbLine in dbLines:
    #    fs.write(dbLine[0][:-1] + "\n")
    fs.write("\n")

    fs.write(">taxonSampling_file\n")
    if os.path.isfile(eachDirAddress + "000_uploadedTxSampling.txt"):
        fs.write("Uploaded\n")
    else:
        fs.write("Not uploaded\n")
    fs.write("\n")


    fs.write(">E-value_threshold_for_reported_sequences\n"  + str(blastEvalue)  + "\n\n")
    fs.write(">Number_of_hits_to_report_per_genome\n" + str(blastTopHits) + "\n\n")
    fs.write(">TreeSearchMethod\n"   + "Neighbor-joining method (Saitou and Nei 1986)\n\n")

    if querySeqChr == "nucl":
        if dataset == "Exclude3rd":
            fs.write(">SubstitutionModel\n" + "TN93 (Tamura and Nei 1993) + gamma\n\n")
        elif dataset == "Include3rd":
            fs.write(">SubstitutionModel\n" + "TN93 (Tamura and Nei 1993) + gamma\n\n")
        else:
            fs.write(">SubstitutionModel\n" + "WAG (Whelan and Goldman 2001) + gamma\n\n")
    else:
        fs.write(">SubstitutionModel\n" + "WAG (Whelan and Goldman 2001) + gamma\n\n")

    fs.write(">Aligned_site_rate_threshold_in_Unambiguously_aligned_sites\n" + str(nonGapSiteRate) + "\n\n")

    fs.write(">Dataset\n" + dataset +  "\n")
    fs.write("\n")

    if reconciliation == "treeSearchRearrange":
        fs.write(">Rearrangement_BS_value_threshold\n")
        fs.write(str(RearrangementBSthreshold) + "\n")
        fs.write("\n")

    fs.write(">Dependencies\n")
    fs.write("ORTHOSCOPEv1.5.2\n")
    fs.write("BLAST 2.7.1+\n")
    fs.write("MAFFT v7.356b\n")
    fs.write("PAL2NAL v13\n")
    fs.write("trimAl 1.2rev59\n")
    fs.write("ape in R, Version: 5.6.2\n")
    if querySeqChr == "prot":
        fs.write("FastME 2.0\n")
    if reconciliation == "treeSearchRearrange":
        fs.write("Notung-2.9\n")
    fs.write("\n")


    fs.close()
    
### makeSummary End


def uniqueList(listAll):
    names_uniq = []
    for x in listAll:
        if x not in names_uniq:
            names_uniq.append(x)
    return names_uniq


def select_blastHitUnique(blastResFileFN):
    f = open(eachDirAddress + blastResFileFN)
    blastResAllLines = list(f)
    f.close()
    nameLinesAll = []
    nameLineAdd = ""
    flag = 0
    for line in blastResAllLines:
        line = line.rstrip("\n")
        if line.startswith("Length"):
            nameLineAdd = re.sub("^> ", ">", nameLineAdd)
            nameLinesAll.append(nameLineAdd)
            nameLineAdd = ""
            flag = 0
        if flag == 1:
            nameLineAdd += line
        if line.startswith(">"):
            flag = 1
            nameLineAdd += line
    return uniqueList(nameLinesAll)


def hitRecPicker():
    blastResOut = open(eachDirAddress + "010_blastRes.txt", "w")
    AAout = open(eachDirAddress + "030_retrievedAAfas.txt", "w")
    CDNAout = ""
    if querySeqChr == "nucl":
        CDNAout = open(eachDirAddress + "030_retrievedDNAfas.txt", "w")

    #'''
    num_BlastpHits = 0
    for dbline in dbLines:

        recdbAA = readFasta_dict(dbAddress, dbline[1])
        recdbDNA = ""
        if querySeqChr == "nucl":
            recdbDNA = readFasta_dict(dbAddress, dbline[2])

        #blastResFile   = dbline[3]
        blastResFile   = "005_vs" + dbline[0][:-1] + ".txt"
        #print("blastResFile:", blastResFile, "<br>")
        uniqueHitLines = select_blastHitUnique(blastResFile)

                
        for line in reversed(uniqueHitLines):
            if not line.startswith(">"):
                continue
            ################################### Blast hit delete START
            #if re.search(">LanceletC_BRBE143570F-t", line):
            #    #print("LanceletC_BRBE143570F-t skipped.")
            #    continue
            #if re.search(">LanceletC_BRBE239040R-t1", line):
            #    continue
            #if re.search("> Drosophila_FBpp0309625_Retinal-Homeobox", line):
            #    continue
            ################################### Blast hit delete END
            
            num_BlastpHits += 1

            line = line.rstrip("\n")
            lineTMP = re.sub(">", ">" + dbline[0], line)
            blastResOut.write(lineTMP + "\n")
            match = re.search("DBNLINE\|(\d+)\|", line)
            DBNLINEnum = match.group(1)
            #print "DBNLINEnum:", DBNLINEnum, "<br>"
            #print list(recdbAA.keys())[int(DBNLINEnum)], "<br><br>"
            name_AA = list(recdbAA.keys())[int(DBNLINEnum)]
            name_AA = re.sub(">", ">" + dbline[0][:-1] + "_", name_AA)
            
            #print("name_AA", name_AA, "<br>")
            name_AA = re.sub("(>.{" + lengthLimit_nameLine + "}).*", r"\1", name_AA)
            #print("name_AA", name_AA, "<br><br>")

            AAout.write(name_AA + "\n")
            AAout.write(list(recdbAA.values())[int(DBNLINEnum)]  + "\n")
            if querySeqChr == "nucl":
                nameLine_CDNA = list(recdbDNA.keys())[int(DBNLINEnum)]
                nameLine_CDNA = re.sub(">", ">" + dbline[0][:-1] + "_", nameLine_CDNA)
                #print("nameLine_CDNA",nameLine_CDNA)
                #exit()
                #CDNAout.write(list(recdbDNA.keys())[int(DBNLINEnum)]   + "\n")
                CDNAout.write(nameLine_CDNA + "\n")
                CDNAout.write(list(recdbDNA.values())[int(DBNLINEnum)] + "\n")
            
            #print()
    blastResOut.close()

    #print "num_BlastpHits:", num_BlastpHits
    if num_BlastpHits < 4:
        result = "Blastp hits was less than 4."
        error_makeSummary(result)
        error_resHtmlMaker(result)
        htmlAddress = make_htmlAddress()
        print ("Blastp hits was less than 4. Cannot estimate gene tree.")
        print (htmlAddress)
        exit()

    recsAA_queryFile   = readFasta_dict(eachDirAddress, "000_aaQuery.txt")
    recsCDNA_queryFile = ""
    if querySeqChr == "nucl":
        recsCDNA_queryFile = readFasta_dict(eachDirAddress, "000_cDNAQuery.txt")

    if reconciliation == "treeSearchRearrange":

        spNames_IndbLines = make_spNames_IndbLines()
        for i in reversed(range(len(recsAA_queryFile))):
            nameLine = list(recsAA_queryFile.keys())[i]
            #print("nameLine:", nameLine)
            match = re.search(">([^_]+)_([^_]+)_", nameLine)
            spName = match.group(1)
            #protID = match.group(2)
            if not spName in spNames_IndbLines:
                #print("spNamePrefix_queryFile:", spNamePrefix_queryFile)
                #print("protID_queryFile:", protID_queryFile)
                #AAout.write(">" + spName + "_" + protID + "_NONE cds gene:" + protID + " transcript:" + protID + " \n")
                AAout.write(nameLine + " \n")
                AAout.write(list(recsAA_queryFile.values())[i] + "\n")
                #CDNAout.write(">" + protID + " cds gene:" + protID + " \n")
                if querySeqChr == "nucl":
                    CDNAout.write(nameLine + " \n")
                    CDNAout.write(list(recsCDNA_queryFile.values())[i] + "\n")
    else:
        for i in reversed(range(len(recsAA_queryFile))):
            #nameLine = list(recsAA_queryFile.keys())[i]
            #match = re.search(">([^_]+)_([^_]+)_", nameLine)
            #spName = match.group(1)
            #protID       = match.group(2)
            #AAout.write(">" + spName + "_" + protID + "_NONE cds gene:" + protID + " transcript:" + protID + " \n")
            AAout.write(list(recsAA_queryFile.keys())[i] + "\n")
            AAout.write(list(recsAA_queryFile.values())[i] + "\n")
            #CDNAout.write(">" + protID + " cds gene:" + protID + " \n")
            if querySeqChr == "nucl":
                CDNAout.write(list(recsCDNA_queryFile.keys())[i] + "\n")
                CDNAout.write(list(recsCDNA_queryFile.values())[i] + "\n")

    AAout.close()
    if querySeqChr == "nucl":
        CDNAout.close()


def make_spNames_IndbLines():
    spNames_IndbLines = []
    for dbLine in dbLines:
        spNamePrefix_dbLine = dbLine[0]
        spNamePrefix_dbLine = re.sub("_$", "", spNamePrefix_dbLine)
        spNames_IndbLines.append(spNamePrefix_dbLine)
    return spNames_IndbLines


def blastpSearch():
    DirAafileName = eachDirAddress + "000_aaQuery.txt"

    for dbline in dbLines:
        dbAAfile = dbAddress + dbline[1]
        #outFile  = eachDirAddress + dbline[3]
        outFile = eachDirAddress + "005_vs" + dbline[0][:-1] + ".txt"
        comLine = "orthoscopeScripts/blastp -query {0} -evalue {1} -num_alignments {2} -num_descriptions {3} -db {4} -out {5}"\
                  .format(DirAafileName, blastEvalue, blastTopHits, blastTopHits, dbAAfile, outFile)

        #print("comLine:", comLine)
        #exit()
        subprocess.call(comLine, shell=True)


def blastResQueryHitPicker(nameLine):
    #print("nameLine2:", nameLine)
    match = re.search(">([^_]+)_", nameLine)
    speciesNamePrefix = match.group(1)

    topHit = ""
    identity = ""

    spNames_IndbLines = make_spNames_IndbLines()

    if not speciesNamePrefix in spNames_IndbLines:
        topHit = nameLine
    else:
        blastOutFileName = "005_vs" + speciesNamePrefix + ".txt"
        #print("blastOutFileName:", blastOutFileName)
        f = open(eachDirAddress + blastOutFileName)
        lines = list(f)
        f.close()
        flag = 0
        for line in lines:

            spNameProtID = nameLine
            #print spNameProtID, "<br>"
            spNameProtID = re.sub(">([^_]+_[^_]+)_.*$", r"\1", spNameProtID)
            key = "^Query= " + spNameProtID

            if flag == 1:
                if line.startswith(">"):
                    topHit = line
                    topHit = re.sub("^> ",">" + speciesNamePrefix + "_", topHit)
                    topHit = re.sub("^> ",">", topHit)
                    topHit = re.sub(" .*$","", topHit)
                    topHit = topHit.rstrip("\n")

                if line.startswith(" Identities ="):
                    match = re.search("(Identities =[^,]+),", line)
                    identity = match.group(1)
                    break

            if re.search(key, line):
                flag = 1

        if not topHit:
            topHit = ">" + speciesNamePrefix + "_XXXX_NoBlastHit"
    #print("topHit", topHit, "<br>")
    topHit_mod = re.sub("(.{" + lengthLimit_nameLine + "}).*", r"\1", topHit)
    #print("topHit_mod", topHit_mod, "<br><br>")
    return topHit_mod, identity


def make_querySeqLines4summary(uploadedSEQfn):
    topHits = topHitPicker()
    topHitValue0s      = [x[0] for x in topHits.values()]
    longestNameInDBLen = len(max(topHitValue0s, key = len))
    #longestuploadedCDNAfnLen = len(max(uploadedSEQfn.keys(), key = len))
    #print("uploadedSEQfn",uploadedSEQfn)
    #print("len(topHits)", len(topHits))
    linesFN = []
    for i in range(len(topHits)):
        uploadedSeqName = list(uploadedSEQfn.keys())[i][1:]
        #print("uploadedSeqName:", uploadedSeqName)
        #print("list(topHits.values())[i][1:]:", list(topHits.values())[i][1:])
        #print("\n")
        uploadedSeqName = re.sub("[\n\r]", "", uploadedSeqName)
        #print("list(topHits.values())[i][0]", list(topHits.values())[i][0])
        #exit()
        if reconciliation == "treeSearchRearrange":
            InfoIdentity = list(topHits.values())[i][1]
            #print("InfoIdentity",InfoIdentity)
            if not InfoIdentity:
                InfoIdentity = "Name"
            #fs.write(list(topHits.values())[i][0][1:] \
            linesFN.append(list(topHits.values())[i][0][1:] \
                     + " " * (longestNameInDBLen - len(list(topHits.values())[i][0][1:])) \
                     + " <= " \
                     + "[" + InfoIdentity + "] "\
                     + " " * (30 - len(InfoIdentity)) \
                     + ": "+ uploadedSeqName \
                     #+ " " * (longestuploadedCDNAfnLen - len(uploadedSeqName)) \
                     + "\n")
        else:
            InfoIdentity = "Name"
            #fs.write(list(topHits.values())[i][0][1:] \
            linesFN.append(list(topHits.values())[i][0][1:] \
                     + " " * (longestNameInDBLen - len(list(topHits.values())[i][0][1:])) \
                     + " <= " \
                     + "[" + InfoIdentity + "] "\
                     + ": "+ uploadedSeqName \
                     + "\n")
    return linesFN


def make_line_for_topHits(uploadedSEQfn, topHits):
    line = ""

    topHitValue0s = [x[0] for x in topHits.values()]
    longestNameInDBLen = len(max(topHitValue0s, key = len))

    for i in range(len(topHits)):
        uploadedSeqName = list(uploadedSEQfn.keys())[i][1:]
        #print("uploadedSeqName:", uploadedSeqName)
        #print("list(topHits.values())[i][1:]:", list(topHits.values())[i][1:])
        #print("\n")
        uploadedSeqName = re.sub("[\n\r]", "", uploadedSeqName)
        #print("list(topHits.values())[i][0]", list(topHits.values())[i][0])
        #exit()

        if reconciliation == "treeSearchRearrange":
            InfoIdentity = list(topHits.values())[i][1]
            if not InfoIdentity:
                InfoIdentity = "Name"
            line += list(topHits.values())[i][0][1:] \
                     + " " * (longestNameInDBLen - len(list(topHits.values())[i][0][1:])) \
                     + " <= " \
                     + "[" + InfoIdentity + "] "\
                     + " " * (30 - len(InfoIdentity)) \
                     + ": "+ uploadedSeqName \
                     + "\n"
        else:
            InfoIdentity = "Name"
            line += list(topHits.values())[i][0][1:] \
                     + " " * (longestNameInDBLen - len(list(topHits.values())[i][0][1:])) \
                     + " <= " \
                     + "[" + InfoIdentity + "] "\
                     + ": "+ uploadedSeqName \
                     + "\n"

    return line



def make_line_for_topHits_html(uploadedSEQfn, topHits):
    line = ""

    topHitValue0s = [x[0] for x in topHits.values()]
    longestNameInDBLen = len(max(topHitValue0s, key = len))

    for i in range(len(topHits)):
        uploadedSeqName = list(uploadedSEQfn.keys())[i][1:]
        #print("uploadedSeqName:", uploadedSeqName)
        #print("list(topHits.values())[i][1:]:", list(topHits.values())[i][1:])
        #print("\n")
        uploadedSeqName = re.sub("[\n\r]", "", uploadedSeqName)
        #print("list(topHits.values())[i][0]", list(topHits.values())[i][0])
        #exit()

        if reconciliation == "treeSearchRearrange":
            InfoIdentity = list(topHits.values())[i][1]
            if not InfoIdentity:
                InfoIdentity = "Name"
            line += list(topHits.values())[i][0][1:] \
                     + " " * (longestNameInDBLen - len(list(topHits.values())[i][0][1:])) \
                     + "<br>&emsp;<= " \
                     + "[" + InfoIdentity + "] "\
                     + " " * (30 - len(InfoIdentity)) \
                     + ": "+ uploadedSeqName \
                     + "<br>"
        else:
            InfoIdentity = "Name"
            line += list(topHits.values())[i][0][1:] \
                     + " " * (longestNameInDBLen - len(list(topHits.values())[i][0][1:])) \
                     + "<br>&emsp;<= " \
                     + "[" + InfoIdentity + "] "\
                     + ": "+ uploadedSeqName \
                     + "<br><br>"
    return line


def topHitPicker():
    topHits = OrderedDict()

    recs_SEQ = ""
    if querySeqChr == "nucl":
        recs_SEQ = readFasta_dict(eachDirAddress, "000_cDNAQuery.txt")
    else:
        recs_SEQ = readFasta_dict(eachDirAddress, "000_aaQuery.txt")
    for nameLine, seq in recs_SEQ.items():
        if reconciliation == "treeSearchRearrange":
            #print "nameLine:", nameLine, "<br>"
            blastTopHitNameLine, blastTopHitIdentity = blastResQueryHitPicker(nameLine)
            #print "blastTopHitIdentity:", blastTopHitIdentity, "<br>"
            topHits[nameLine] = [blastTopHitNameLine, blastTopHitIdentity.rstrip("\n")]
        else:
            nameLineTMP = re.sub(" .*$","", nameLine)
            #topHits[nameLine] = nameLineTMP
            topHits[nameLine] = [nameLineTMP, 1]
    return topHits


### Blast End


### seqSelect_nonGapSitesEnoughLong Start
def compare2seqs(querySeq, otherSeq):
    count = 0
    for i in range(len(querySeq)):
        #if re.match(querySeq[i], "[-Xx]") and re.match(querySeq[i], "[-Xx]"):
        #    print("  NOT")
        #elif re.match(querySeq[i], "[-Xx]") and re.match(querySeq[i], "\w"):
        #    print("  COUNT")
        #    count += 1
        #if re.match(querySeq[i], "\w") and re.match(querySeq[i], "\w"):
        if re.match("\w", querySeq[i]) and re.match("\w", otherSeq[i]):
            #print("  COUNT")
            count += 1
        elif re.match("-", querySeq[i]) and re.match("\w", otherSeq[i]):
            #print("querySeq[i]:", querySeq[i])
            #print("otherSeq[i]:", otherSeq[i])
            #print("  COUNT")
            count += 1
        else:
            #print("querySeq[i]:", querySeq[i])
            #print("otherSeq[i]:", otherSeq[i])
            #print("  NOT")
            continue
    return round(float(count)/len(querySeq),3)


def calculate_nonGapSiteRate(trimledAAfile):
    recs_nonGapSiteRateFN = OrderedDict()
    recAA = readFasta_dict(eachDirAddress, trimledAAfile)

    querySeq = list(recAA.values())[-1]
    for otherName, otherSeq in recAA.items():
        protID = re.sub(" .*$","", otherName)
        aaRate = compare2seqs(querySeq, otherSeq)
        recs_nonGapSiteRateFN[protID] = aaRate
    return recs_nonGapSiteRateFN


def seqSelect_nonGapSitesEnoughLong():

    recs_nonGapSiteRate = calculate_nonGapSiteRate("042_AA.fas.trm")

    file_overRateAA = open(eachDirAddress + "044_overRateAA.fas", "w")
    file_overRateDNA = ""
    if querySeqChr == "nucl":
        file_overRateDNA = open(eachDirAddress + "044_overRateDNA.fas", "w")
    file_unambSiteRatefile = open(eachDirAddress + "044_unambSiteRate.txt", "w")

    recAA = readFasta_dict(eachDirAddress, "030_retrievedAAfas.txt")
    recDNA = ""
    if querySeqChr == "nucl":
        recDNA = readFasta_dict(eachDirAddress, "030_retrievedDNAfas.txt")
    #print("recs_nonGapSiteRate", len(recs_nonGapSiteRate))
    #print("recAA.keys()", str(len(recAA)))
    #print("recDNA.keys()", str(len(recDNA)))
    #exit()
    for i in range (len(recs_nonGapSiteRate)):
        name = list(recs_nonGapSiteRate.keys())[i]
        rate = list(recs_nonGapSiteRate.values())[i]
        #print(i, rate, "<br>")
        if rate > nonGapSiteRate:
            file_overRateAA.write  (list(recAA.keys())[i]    + "\n")
            file_overRateAA.write  (list(recAA.values())[i]  + "\n")
            if querySeqChr == "nucl":
                nameLineTMP = list(recDNA.keys())[i]
                #nameLineTMP = re.sub(" +.*$", "", nameLineTMP)
                #print i, nameLineTMP, "<br>"
                file_overRateDNA.write(list(recDNA.keys())[i]   + "\n")
                file_overRateDNA.write(list(recDNA.values())[i] + "\n")
        file_unambSiteRatefile.write(name + "\n" + str(rate) + "\n")
    file_overRateAA.close()
    if querySeqChr == "nucl":
        file_overRateDNA.close()
    file_unambSiteRatefile.close()

    
### seqSelect_nonGapSitesEnoughLong End


### alignmentFile_PhyAnal Start
def orderedDict2FasFile(recs, outfile):
    out           = open(eachDirAddress + outfile, "w")
    for name,value in recs.items():
        out.write(name + "\n")
        out.write(value + "\n")
    out.close()

def orderedDict2phyFile(recs, outfile):
    secLength     = len(sorted(recs.values())[0])
    spSeqSizeLine = str(len(recs)) + " " + str(secLength)

    recs = whiteSpaceAdd(recs)
    out           = open(eachDirAddress + outfile, "w")
    out.write(spSeqSizeLine + "\n")
    for name,value in recs.items():
        out.write(name + value + "\n")
    out.close()


def read_TrimalHTMLout_dict(trimAAresult):
    f = open(eachDirAddress + trimAAresult)
    lines = list(f)
    f.close()

    recs_trimalHTMLout  = OrderedDict()
    for line in lines:
        if line.startswith("    <span class=sel>"):
            line = line.rstrip("\n")
            match = re.search("<span class=sel>([^<]+)</span> +([^ ].*)$", line)
            name     = match.group(1)
            sequence = match.group(2)
            if not name in recs_trimalHTMLout.keys():
                recs_trimalHTMLout[name] = ""
            recs_trimalHTMLout[name] += sequence

    #for name, seq in recs_trimalHTMLout.items():
    #    print(name)
    #    print(seq)
    #exit()
    return recs_trimalHTMLout


def trimaledFileMakerDNA(fastaFile, trimAAresult):
    recs = readFasta_dict(eachDirAddress, fastaFile)
    #for name, seq in recs.items():
    #    print(name)
    #    print(seq)
    
    recs_trimalHTMLout    = read_TrimalHTMLout_dict(trimAAresult)
    trimalMarkedSites = list(recs_trimalHTMLout.values())[-1]
    trimalMarkedSites = re.sub("<span class=sel>.</span>", "#", trimalMarkedSites)

    recsTrimed = OrderedDict()
    for name,value in recs.items():
        sequence = ""
        for i in range(len(trimalMarkedSites)):
            if trimalMarkedSites[i] == "#":
                #out.write(trimalMarkedSites[i])
                sequence += value[i*3] + value[i*3+1] + value[i*3+2]
        recsTrimed[name] = sequence
    orderedDict2phyFile(recsTrimed, outfile = "080_trimedCDNAPhy.txt")

    #recsTrimedExc3rd = OrderedDict()
    #for name,value in recsTrimed.items():
    #    seqeucneExc3rd = ""
    #    for i in range(len(value)):
    #        if i % 3 == 0 or i % 3 == 1:
    #            seqeucneExc3rd += value[i]
    #    recsTrimedExc3rd[name] = seqeucneExc3rd
    #orderedDict2phyFile(recsTrimedExc3rd, outfile = "080_trimeExc3rdPhy.txt")

### alignmentFile_PhyAnal End


### alignmentFile_HTML
def outGroupSelect (phyFileName):
    recSeqFN = readPhy_dict(phyFileName)
    outgroupTMP = list(recSeqFN.keys())[0]
    return outgroupTMP[1:]

def reorderSeqByTree(recsFN, treeFileName):
    leaves = []
    leaves = leafCollectInOrderFrom_bothNHXnewick(treeFileName)
    seqDictFN  = OrderedDict()
    for leaf in reversed(leaves):
        Lleaf = ">" + leaf
        seqDictFN[Lleaf] = recsFN[Lleaf]
    return seqDictFN

def trimAlresReader(inFileAAFN):
    aafile = open(eachDirAddress + inFileAAFN, "r")
    aafileLines = list(aafile)
    aafile.close()
    
    sequenceTMP = ""
    name_firstRec = ""
    for line in aafileLines:
        if re.search("    <span class=sel>[^<]+<", line):
            match = re.search("    <span class=sel>([^<]+)<", line)
            name_firstRec = match.group(1)
            break
    for line in aafileLines:
        if re.search("<span class=sel>" + name_firstRec + "<", line):
            line = re.sub("^ +<span class=sel>" + name_firstRec + "</span> +", "", line).rstrip("\n")
            line = re.sub("<span class=sel>.</span>", "#", line)
            sequenceTMP += line

    sequence = ""
    for chr in sequenceTMP:
        codon = ""
        if chr == "#":
            codon = "111"
        else:
            codon = "000"
        sequence += codon

    if sequence:
        return sequence
    else:
        print("No key sequence line in :", inFileAAFN)
        print("  searched by :", topHitName_1stQuery[1:])
        exit()

def addSiteNotes(inFileAAFN, recsFN):
    seqDictFN  = OrderedDict()
    sequence_tm = trimAlresReader(inFileAAFN)
    seqDictFN[">Used4treeSearch"] = sequence_tm
    for name, value in recsFN.items():
      seqDictFN[name] = value
    return seqDictFN
    
#def speciesColorStyleAdd(recsFN1):
#    recsFN2 = OrderedDict()
#    for name, value in recsFN1.items():
#        nameStyled = ""
#        for dbEle in dbLines:
#            if re.search(dbEle[0], name):
#                nameStyled = "<span class=" + dbEle[5] +">" + name[1:] + "</span>"
#                break
#        if nameStyled:
#            recsFN2[nameStyled] = value
#        else:
#            recsFN2[name[1:]] = value
#    return recsFN2

def delete_nameSpaceSeqBp(recsFN):
    recsFN2        = OrderedDict()
    for name, seq in recsFN.items():
        if re.search(" \d+ bp$", name):
            name = re.sub(" \d+ bp$", "", name)
        recsFN2[name] = seq
    return recsFN2

def arrowheadDelete(recsFN):
    recsFN2        = OrderedDict()
    for name, seq in recsFN.items():
        if name.startswith(">"):
            recsFN2[name[1:]] = seq
        else:
            recsFN2[name] = seq
    return recsFN2

def nameChange_whiteLaterDelete(recsFN):
    recsFN2        = OrderedDict()
    for name, seq in recsFN.items():
        name = re.sub(" .*$", "", name)
        recsFN2[name] = seq
    return recsFN2

def gapDelete(recsFN):
    recsFN2        = OrderedDict()
    for name, seq in recsFN.items():
        seq = re.sub("-", "", seq)
        recsFN2[name] = seq
    return recsFN2


def whiteSpaceAdd(recsFN1):
    longestName = max(recsFN1.keys(), key = len)
    longestName = re.sub("<[^>]+>", "", longestName)
    longestNameLen = len(longestName)
    recsFN2        = OrderedDict()
    for name,value in recsFN1.items():
        name = re.sub("^>", "", name)
        nameTMP = re.sub("<[^>]+>", "", name)
        nameWhiteSpace = name + " " * (longestNameLen - len(nameTMP) + 2)
        recsFN2[nameWhiteSpace] = value
    return recsFN2

def keySequencePicker(topHitName_1stQuery, recsFN):
    hitNames = [name for name in recsFN.keys() if re.search(topHitName_1stQuery[1:], name)]
    if len(hitNames) > 1:
        print("More than 2 keySpecies:%s. Stopped." % topHitName_1stQuery[1:])
        exit()
    elif len(hitNames) == 0:
        print("For your query sequences, no BLAST hit was found for the first query species: %s." % topHitName_1stQuery[1:])
        print("Please try Estimating gene tree mode.")
        exit()
    else:
        return recsFN[hitNames[0]]


def matchFirsterCodon (topHitName_1stQuery, keyCodonFN, nameFN, codonFN):
    codonFNmod = ""
    nameFN = ">" + nameFN
    if re.search(topHitName_1stQuery, nameFN):
        codonFNmod = codonFN
    else:
        for i in range(0,3):
            if keyCodonFN[i] == "X" or keyCodonFN[i] == "-" or keyCodonFN[i] == "N":
                codonFNmod += codonFN[i]
            elif keyCodonFN[i] == codonFN[i]:
                codonFNmod += "."
            else:
                codonFNmod += codonFN[i]
    return codonFNmod


def codonHTMLpicker(topHitName_1stQuery, keySpCodonFN, nameFN, codonFN):
    codonMatchFirsted = ""
    
    if re.search("\d\d\d", codonFN) or re.search("-", codonFN):
        codonMatchFirsted = codonFN
    else:
        codonMatchFirsted = geneticCodeHTML.get(codonFN, "CODON")
        codonFN2          = matchFirsterCodon(topHitName_1stQuery, keySpCodonFN, nameFN, codonFN)
        codonMatchFirsted = re.sub("CODON", codonFN2, codonMatchFirsted)
    return codonMatchFirsted


def recOneLinesMaker(topHitName_1stQuery, keySequenceFN, startPosNF, stopPosNF, recsFN):

    htmlOneLines = []

    longestName    = max(recsFN.keys(), key = len)
    longestName    = re.sub("<[^>]+>", "", longestName)
    longestNameLen = len(longestName)
    NumberLine     = " " * (longestNameLen) + str(startPosNF+1) + "\n"
    htmlOneLines.append(NumberLine)

    for name in recsFN.keys():
        htmlOneLine = name
        for p in range(startPosNF, stopPosNF, 3):
            if p+2 <= stopPosNF:
                keySpCodon   = keySequenceFN[p] + keySequenceFN[p+1] + keySequenceFN[p+2]
                codon        = recsFN[name][p]  + recsFN[name][p+1]  + recsFN[name][p+2]
                codonHTML    = codonHTMLpicker(topHitName_1stQuery, keySpCodon, name, codon)
                htmlOneLine += codonHTML

        spName = re.sub("[^_]+_[^_]+$", "", name)
        htmlOneLine = "<span class=" + spName + ">" + htmlOneLine + "</span>"
        if re.search(topHitName_1stQuery[1:], name):
            htmlOneLine = "<u>" + htmlOneLine + "</u>"

        htmlOneLine += "\n"
        htmlOneLines.append(htmlOneLine)

    return htmlOneLines


def make_alignHTMLtop():
    styles = styles_html_codons
    for spNameTMP in taxonSamplingList:
        match = re.search("^([^_]+_)([^_]+)$", spNameTMP)
        spName = match.group(1)
        color = match.group(2)
        styles += ("                ." + spName + " { color: " + color + ";  }\n")
    alignHTMLtopTMP1 = re.sub("STYLELINES", styles, alignHTMLtopTMP)
    return alignHTMLtopTMP1

def codonHTMLfileMaker(oneLineLengthFN, recsFN, outFileFN):

    topHits = topHitPicker()
    topHitName_1stQuery = list(topHits.values())[0][0]
    topHitSeq_1stQuery  = keySequencePicker(topHitName_1stQuery, recsFN)

    siteNum  = 1
    startPos = 0
    stopPos  = 0
    sequenceLength = len(topHitSeq_1stQuery)
    #print("sequenceLength:", sequenceLength)
    
    alignHTMLtop = make_alignHTMLtop()
    out = open(eachDirAddress + outFileFN, "w")
    out.write(alignHTMLtop)
    out.write(str(len(recsFN)) + " " + str(sequenceLength) + "\n")
    for i in range(sequenceLength):
        if i > 1 and i % oneLineLengthFN == 0:
            startPos = i - oneLineLengthFN
            stopPos = startPos + oneLineLengthFN
            #print("siteNum :", startPos + 1)
            #print("startPos:", startPos)
            #print("stopPos :", stopPos)
            #for j in range(startPos, stopPos):
            #    print(str(j) + ",", end="")
            #print()
            recOneLines = recOneLinesMaker(topHitName_1stQuery, topHitSeq_1stQuery, startPos, stopPos, recsFN)
            for line in recOneLines:
                out.write(line)
            out.write("\n")

    #print("siteNum :", stopPos + 1)
    #print("startPos:", stopPos)
    #print("stopPos :", sequenceLength)
    #for i in range(stopPos, sequenceLength):
    #    print(str(i) + ",", end="")
    #print()
    recOneLines = recOneLinesMaker(topHitName_1stQuery, topHitSeq_1stQuery, stopPos, sequenceLength, recsFN)
    for line in recOneLines:
        out.write(line)
    out.write(alignHTMLbottom)
    out.close()

def cDNAaligHtmlMaker(inFile, inFileAA, tree4leafOrder, outFile):
    oneLineLength = 90         ## For cDNA, needs to be a multiple of 3
    recs        = readFasta_dict(eachDirAddress, inFile)
    recs        = reorderSeqByTree(recs, tree4leafOrder)
    recs        = addSiteNotes(inFileAA, recs)
    recs        = whiteSpaceAdd(recs)
    codonHTMLfileMaker(oneLineLength, recs, outFile)

def reorderDeleteGap_Fas2FasByTree(fastaFileName, tree4leafOrder, outPhyFileName):
    recs = readFasta_dict(eachDirAddress,fastaFileName) 
    recs = nameChange_whiteLaterDelete(recs)
    recs = gapDelete(recs)
    recs = reorderSeqByTree(recs, tree4leafOrder)
    orderedDict2FasFile(recs, outfile = outPhyFileName)

def phy2fastmePhy(phyFileName, outFastmePhyFileName):
    recsFN = readPhy_dict(phyFileName)
    secLength     = len(sorted(recsFN.values())[0])
    spSeqSizeLine = str(len(recsFN)) + " " + str(secLength)

    recsFN = whiteSpaceAdd(recsFN)
    out  = open(eachDirAddress + outFastmePhyFileName, "w")
    out.write(spSeqSizeLine + "\n")
    for name,value in reversed(recsFN.items()):
        out.write(name + value + "\n")
    out.close()


def fas2phy(fastaFileName, outPhyFileName):
    recs = readFasta_dict(eachDirAddress,fastaFileName)
    recs = delete_nameSpaceSeqBp(recs)
    orderedDict2phyFile(recs, outfile = outPhyFileName)

def reorderFas2PhyByTree(fastaFileName, tree4leafOrder, outPhyFileName):
    recs = readFasta_dict(eachDirAddress,fastaFileName) 
    recs = delete_nameSpaceSeqBp(recs)
    for name, seq in recs.items():
        print(name)
        print("<br><br>")
    exit()
    recs = reorderSeqByTree(recs, tree4leafOrder)
    #recs = arrowheadDelete(recs)
    orderedDict2phyFile(recs, outfile = outPhyFileName)

def reorderPhyFileByTree(phyFileName, tree4leafOrder, outPhyFileName):
    recs = readPhy_dict(phyFileName) 
    #print(recs.values())
    recs = reorderSeqByTree(recs, tree4leafOrder)
    #recs = arrowheadDelete(recs)
    #recs = whiteSpaceAdd(recs)
    orderedDict2phyFile(recs, outfile = outPhyFileName)

def codonSepalate(recs):
    recs1 = OrderedDict()
    recs2 = OrderedDict()
    recs3 = OrderedDict()
    for name, sec in recs.items():
        recs1[name] = ""
        recs2[name] = ""
        recs3[name] = ""
        for i in range(len(sec)):
            if   i%3 == 0:
               recs1[name] += sec[i]
            elif i%3 == 1:
               recs2[name] += sec[i]
            else:
               recs3[name] += sec[i]
    return recs1, recs2, recs3

def phyCodonToBlock(phyFileName, blockNum, outfile):
    recs = readPhy_dict(phyFileName)
    recs1, recs2, recs3 = codonSepalate(recs)
    recsS1 = OrderedDict()
    if blockNum == 3:
        for name in recs1.keys():
            recsS1[name] = recs1[name] + recs2[name] + recs3[name]
            
    if blockNum == 2:
        for name in recs1.keys():
            recsS1[name] = recs1[name] + recs2[name]

    #recsS1   = arrowheadDelete(recsS1)
    recsS2 = whiteSpaceAdd(recsS1)

    out = open(eachDirAddress + outfile, "w")
    out.write(str(len(recsS2)) + " " + str(len(list(recsS2.values())[0])) + "\n")
    for name,sec in recsS2.items():
        out.write(name + sec + "\n")
    out.close()


### alignmentFile_HTML End

### result_HTML Start
def error_resHtmlMaker(resultFN):
    topHits = topHitPicker()
    topHitName_1stQuery = list(topHits.values())[0]
    resHTMLlines1 = re.sub('FIRSTQUERY', topHitName_1stQuery[0][1:], resHTMLlines_error)
    resHTMLlines1 = re.sub('RESULTMESSAGE', resultFN, resHTMLlines1)

    if analysisGroupName == "Protostomia":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: protostome orthogroup', resHTMLlines1)
    elif analysisGroupName == "Deuterostomia":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: deuterostome orthogroup', resHTMLlines1)
    elif analysisGroupName == "Acropora":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: acropora orthogroup', resHTMLlines1)
    elif analysisGroupName == "Vertebrata":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: vertebrate orthogroup', resHTMLlines1)
    elif analysisGroupName == "Mammalia":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: Mammalian orthogroup', resHTMLlines1)
    elif analysisGroupName == "Plants":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: Plant orthogroup', resHTMLlines1)
    elif analysisGroupName == "Chordata":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: chordate orthogroup', resHTMLlines1)
    elif analysisGroupName == "Actinopterygii":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: actinopterygian orthogroup', resHTMLlines1)
    else:
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: orthogroup', resHTMLlines1)
    if reconciliation == "treeSearchRearrange":
        resHTMLlines1 = re.sub("REARRANGED_TREE",  '<img src="115_1stRearranged_geneTree.png">', resHTMLlines1)
    else:
        resHTMLlines1 = re.sub(".*REARRANGED.*\n", "", resHTMLlines1)

    out = open(eachDirAddress + "300_resultsREA.html", "w")
    out.write(resHTMLlines1)
    out.close()


def resHtmlMaker():
    
    resHTMLlines1 = ""
    if reconciliation == "treeSearchRearrange":
        f1stTree = open(eachDirAddress + "085_NJBS1st.txt.rearrange.0")
        treeNHX = list(f1stTree)[0]
        f1stTree.close()
        orthoGroup = orthoGroupIdentifyFromNHX (treeNHX)
        
        Res_orthogroupSuppor = make_Res_orthogroupSuppor(orthoGroup[2])
        Res_orthogroupSuppor = 'Result: <br>' +  Res_orthogroupSuppor
        #print("downloadLine:", downloadLine)
        resHTMLlines1 = re.sub('RESULTMESSAGE', Res_orthogroupSuppor, resHTMLlines)
    else:
        resHTMLlines1 = re.sub('.*RESULTMESSAGE.*', '', resHTMLlines)

    downloadFileName = "result" + str(dirName_count) + ".zip"
    downloadMessage = 'Download: <a href="' + downloadFileName + '">' + downloadFileName + '</a>.'
    if reconciliation == "treeSearchRearrange":
        uploadedSEQfn = readFasta_dict(eachDirAddress, "000_uploaded.txt")
        topHits = topHitPicker()
        downloadMessage += "<br><br>Your queries were replaced with the blast top hits:<br>" + make_line_for_topHits_html(uploadedSEQfn, topHits) + "<br>"
    resHTMLlines1 = re.sub('DOWNLOADMESSAGE', downloadMessage, resHTMLlines1)

    if analysisGroupName == "Protostomia":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: protostome orthogroup', resHTMLlines1)
    elif analysisGroupName == "Deuterostomia":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: deuterostome orthogroup', resHTMLlines1)
    elif analysisGroupName == "Acropora":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: acropora orthogroup', resHTMLlines1)
    elif analysisGroupName == "Vertebrata":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: vertebrate orthogroup', resHTMLlines1)
    elif analysisGroupName == "Mammalia":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: mammalian orthogroup', resHTMLlines1)
    elif analysisGroupName == "Plants":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: plant orthogroup', resHTMLlines1)
    elif analysisGroupName == "Chordata":
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: chordate orthogroup', resHTMLlines1)
    else:
        resHTMLlines1 = re.sub("TITLE",  'ORTHOSCOPE: actinopterygian orthogroup', resHTMLlines1)
    if reconciliation == "treeSearchRearrange":
        resHTMLlines1 = re.sub("REARRANGED_TREE",  '<img src="115_1stRearranged_geneTree.png">', resHTMLlines1)
    else:
        resHTMLlines1 = re.sub(".*REARRANGED.*\n", "", resHTMLlines1)
    out = open(eachDirAddress + "300_resultsREA.html", "w")
    out.write(resHTMLlines1)
    out.close()


def make_htmlAddress():
    htmlAddress = ""
    #print("EMAIL_SUBJECT_PREFX", EMAIL_SUBJECT_PREFX)
    if "fish-evol" in EMAIL_SUBJECT_PREFX:
        htmlAddress = '<br>Finished: <a href="http://fish-evol.unit.oist.jp/orthoscopeWork/' + str(dirname_rand) + '/300_resultsREA.html" target="_blank">result' + str(dirname_rand) + '</a><br><br>'
    elif "sakura" in EMAIL_SUBJECT_PREFX:
        htmlAddress = '<br>Finished: <a href="http://153.126.167.45/orthoscopeWork/' + str(dirname_rand) + '/300_resultsREA.html" target="_blank">result' + str(dirname_rand) + '</a><br><br>'
    elif "yamasati" in EMAIL_SUBJECT_PREFX:
        htmlAddress = '<br>Finished: <a href="http://yamasati.nig.ac.jp/orthoscopeWork/' + str(dirname_rand) + '/300_resultsREA.html" target="_blank">result' + str(dirname_rand) + '</a><br><br>'
    elif "yurai.aori.u-tokyo.ac.jp" in EMAIL_SUBJECT_PREFX:
        htmlAddress = '<br>Finished: <a href="http://yurai.aori.u-tokyo.ac.jp/orthoscopeWork/' + str(dirname_rand) + '/300_resultsREA.html" target="_blank">result' + str(dirname_rand) + '</a><br><br>'
    elif "rx1000.site" in EMAIL_SUBJECT_PREFX:
        htmlAddress = '<br>Finished: <a href="http://orthoscope.jp/orthoscopeWork/' + str(dirname_rand) + '/300_resultsREA.html" target="_blank">result' + str(dirname_rand) + '</a><br><br>'
    else:
        print("Check make_htmlAddress.")
        exit()
    return htmlAddress


def compression():
    # in cgi-bin directory
    # by tarfile in cgi-bin directory
    #resDirName = "results" + str(dirname_rand)
    #os.mkdir   (resDirName)
    #shutil.copy(eachDirAddress + "050_retAAseqs.maf",             resDirName + "/aa.txt")
    #shutil.copy(eachDirAddress + "054_p2nOutcDNAfas.txt",                  resDirName + "/cDNA.txt")
    #shutil.copy(eachDirAddress + "087_inroot100Tree.outfile.nwk", resDirName + "/tree.txt")
    #shutil.copy(eachDirAddress + "087_inroot100Tree.outfile.pdf", resDirName + "/tree.pdf")

    # in html directory
    # by make_archive
    resDirName = eachDirAddress + "result" + str(dirName_count)
    os.mkdir   (resDirName)
    shutil.copy(eachDirAddress + "100_1stAnalysisSummary.txt",              resDirName + "/000_summary.txt")
    shutil.copy(eachDirAddress + "100_1stAnalysisSummary.html",              resDirName + "/000_summary.html")
    shutil.copy(eachDirAddress + "042_AA.fas.trm.html",                     resDirName + "/010_aln_prot.html")
    #shutil.copy(eachDirAddress + "087_inroot100Tree.outfile.nwk",          resDirName + "/020_tree.txt")
    if os.path.isfile(eachDirAddress + "115_1stGeneTree.pdf"):
        shutil.copy(eachDirAddress + "115_1stGeneTree.pdf",                     resDirName + "/020_tree.pdf")
    else:
        shutil.copy(eachDirAddress + "115_1stGeneTree.png",                     resDirName + "/020_tree.png")
    if reconciliation == "treeSearchRearrange":
        shutil.copy(eachDirAddress + "115_1stRearranged_geneTree.pdf",      resDirName + "/020_treeRearranged.pdf")

    shutil.copy(eachDirAddress + "240_reorderedAAfas.txt",                   resDirName + "/010_candidates_prot.txt")
    if querySeqChr == "nucl":
        shutil.copy(eachDirAddress + "230_CDNA.html",                        resDirName + "/010_aln_nucl.html")
        shutil.copy(eachDirAddress + "240_reorderedCDNAfas.txt",                 resDirName + "/010_candidates_nucl.txt")

    shutil.copy(dirAddress + "/" + version_orthoscope + "/examples/100_2ndTree.zip",           resDirName + "/100_2ndTree.zip")

    # by zip with glob in html directory
    zipfiles = glob.glob(resDirName + '/*')
    fzip = zipfile.ZipFile(resDirName + '.zip', 'w', zipfile.ZIP_DEFLATED)
    for file in zipfiles:
        fzip.write(file, os.path.basename(file))
    fzip.close()
### result_HTML End


### Count

def makeCount():
    import fcntl
    
    dat = "./orthoscopeScripts/count.dat"
    
    fh = open(dat, "r+")
    fcntl.flock(fh.fileno(), fcntl.LOCK_EX)

    count = fh.read()
    count = int(count)
    count += 1
    count = str(count)

    fh.seek(0)
    fh.write(count)
    fcntl.flock(fh.fileno(), fcntl.LOCK_UN)

    fh.close()
    return count

###


######################################################################################################################
#################################### Main program ####################################################################
######################################################################################################################



#########################
### from HTML
startTime = time.time()
dirName_count = makeCount()
letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
letters_rand = ''.join([random.choice(letters) for _ in range(6)])

dirname_rand = dirName_count + "-" + letters_rand
#dirname_rand = "9476-KlPkhh"
eachDirAddress = dirAddress + "orthoscopeWork/" + str(dirname_rand) + "/"

delete_dirs()

taxonSamplingList, querySeqChr, blastEvalue, blastTopHits, nonGapSiteRate, dataset, RearrangementBSthreshold, analysisGroupName, reconciliation = uploadedFileSave(cgi.FieldStorage())
#print("taxonSamplingList", taxonSamplingList, "\n\n")
#print("\n")
#exit()

#####
### from terminal

'''##### time fixed
#dirname_rand = "311-pVqxob"
##eachDirAddress = "../Documents/orthoscopeWork/" + str(dirname_rand) + "/"
eachDirAddress = "../html/orthoscopeWork/" + str(dirname_rand) + "/"
blastEvalue    = 1e-3
blastTopHits   = 5
nonGapSiteRate = 0.45
#dataset        = "Include3rd"
dataset        = "Exclude3rd"
RearrangementBSthreshold = 60
#taxonSamplingList  = ['Celegans_', 'Silkworm_', 'Drosophila_', 'AcornWormS_', 'AcornWormP_', 'SeaUrchinS_', 'AcanthasterGbr_', 'LanceletF_', 'LanceletC_', 'Oikopleura_', 'Botryllus_', 'CionaS_', 'CionaI_', 'Gar_', 'Medaka_', 'Xenopus_', 'Anole_', 'Chicken_', 'Human_']
taxonSamplingList  = ['Crassostrea-gigas_Black', 'Pinctada-fucata_Black', 'Caenorhabditis-elegans_Black', 'Drosophila-melanogaster_Black', 'Saccoglossus-kowalevskii_Green', 'Ptychodera-flava_Green', 'Acanthaster-planci_Purple', 'Strongylocentrotus-purpuratus_Purple', 'Branchiostoma-belcheri_Orange', 'Branchiostoma-floridae_Orange', 'Oikopleura-dioica_Magenta', 'Ciona-savignyi_Magenta', 'Ciona-intestinalis_Magenta', 'Botryllus-schlosseri_Magenta', 'Gallus-gallus-E102_Blue', 'Homo-sapiens-E102_Blue']
analysisGroupName = "Deuterostomia"
querySeqChr = "nucl"
#reconciliation = "treeSearchOnly"
reconciliation = "treeSearchRearrange"
#####
'''#########################


dbLines = dbLinesMake()
#for dbLine in dbLines:
#    print("dbLine", dbLine, "<br>")
#exit()
nodeNames_ancestral = get_nodeNames_ancestral()

#for dbLine in dbLines:
#  print(dbLine[0])
#  print("<br><br>")
#exit()


'''
'''

check_feasibility_of_completion()
checkUplodedFileAsFastaForamt()

if querySeqChr == "nucl":
    cDNAqueryFileMaker()


allNodes = []
if reconciliation == "treeSearchRearrange":
    allNodes = TEST_querySpIsInTree()

nums_querySeqs = aaSeqMaker()




##print("### blast: 010_blastRes.txt 030_retrievedAAfas.txt 030_retrievedDNAfas.txt<br>")
blastpSearch()
#exit()

hitRecPicker()


##print("### mafft: 040_mafOutAA.txt<br>")
maffLine1 = "orthoscopeScripts/mafft --anysymbol " + eachDirAddress + "030_retrievedAAfas.txt > " + eachDirAddress + "040_mafOutAA.txt"
#print("maffLine1", maffLine1)
subprocess.call(maffLine1, shell=True)


mafftOutFile = open(eachDirAddress + "040_mafOutAA.txt")
mafftOutFileCont = list(mafftOutFile)
if not mafftOutFileCont:
    print ("Error in the amino acid alignment using MAFFT. <br>")
    print ("Probably our database contains U in an amino acid sequence of selected species. <br>Contact jun.inoue@nig.ac.jp.<br>")
    exit()


##print("### trimal: 042_AA.fas.trm<br>")
trimLine1 = "orthoscopeScripts/trimal -out " + eachDirAddress + "042_AA.fas.trm -htmlout " + eachDirAddress + "042_AA.fas.trm.html -in " + eachDirAddress + "040_mafOutAA.txt -gappyout"
#print("trimLine1:", trimLine1, "\n");
subprocess.call(trimLine1, shell=True)
#exit()


##print("### seqSelect_nonGapSitesEnoughLong: 044_overRateAA.fas 044_unambSiteRate.txt 044_overRateDNA.fas<br>")
seqSelect_nonGapSitesEnoughLong()


##print("### mafft: 050_retAAseqs.maf\n")
maffLine2 = "orthoscopeScripts/mafft --anysymbol " + eachDirAddress + "044_overRateAA.fas > " + eachDirAddress + "050_retAAseqs.maf"
##print(maffLine2)
subprocess.call(maffLine2, shell=True)
##print("<br><br>")

#print("### trimal: 052_AA.fas.trm.html<br>")
trimLine2 = "orthoscopeScripts/trimal -out " + eachDirAddress + "052_AA.fas.trm -htmlout " + eachDirAddress + "052_AA.fas.trm.html -in " + eachDirAddress + "050_retAAseqs.maf -gappyout"
#print(trimLine2, "\n");
subprocess.call(trimLine2, shell=True)
##print("<br><br>")



##print("### pal2nal: 054_p2nOutcDNAfas.txt<br>")
if querySeqChr == "nucl":
    pal2nalLine = "orthoscopeScripts/pal2nal.pl " + eachDirAddress + "050_retAAseqs.maf " + eachDirAddress + "044_overRateDNA.fas -output fasta > " + eachDirAddress +"054_p2nOutcDNAfas.txt"
    #print (pal2nalLine)
    subprocess.call(pal2nalLine, shell=True)
    out_p2n = open(eachDirAddress + "054_p2nOutcDNAfas.txt")
    content_out_p2n = list(out_p2n)
    if not content_out_p2n:
        result = "Error in the ORTHOSCOPE database: <br>An inconsistent translation was found. Please report this error to jinoue@g.ecc.u-tokyo.ac.jp"
        print (result + "<br>")
        exit()


if querySeqChr == "nucl":
    trimaledFileMakerDNA("054_p2nOutcDNAfas.txt", "052_AA.fas.trm.html")

    #print("### phyCodonToBlock: 082_trimedBlockInc3rdPhy.txt")
    phyCodonToBlock("080_trimedCDNAPhy.txt", 2, outfile="082_trimedBlockExc3rdPhy.txt")
    phyCodonToBlock("080_trimedCDNAPhy.txt", 3, outfile="082_trimedBlockInc3rdPhy.txt")

    phy2fastmePhy(phyFileName="082_trimedBlockExc3rdPhy.txt", outFastmePhyFileName="082_trimedBlockExc3rdFastmePhy.txt")
    phy2fastmePhy(phyFileName="082_trimedBlockInc3rdPhy.txt", outFastmePhyFileName="082_trimedBlockInc3rdFastmePhy.txt")

fas2phy(fastaFileName="052_AA.fas.trm", outPhyFileName="080_trimedAAPhy.txt")
phy2fastmePhy(phyFileName="080_trimedAAPhy.txt", outFastmePhyFileName="082_trimedAAFastmePhy.txt")



### NJ
outgroup1 = outGroupSelect("080_trimedAAPhy.txt")
if querySeqChr == "nucl":
    if dataset == "Exclude3rd":
        #NJBSline1 = "orthoscopeScripts/fastme -i " + eachDirAddress + "082_trimedBlockExc3rdFastmePhy.txt -d F84 -m BioNJ -b 100 -v 3 -o " + eachDirAddress + "085_NJBS1st.txt > " + eachDirAddress + "085_fastmelog.txt"
        NJBSline1 = "orthoscopeScripts/Rscript orthoscopeScripts/085_NJBSa.R " + eachDirAddress + "082_trimedBlockExc3rdPhy.txt " + outgroup1 + " TN93 " + eachDirAddress + "085_NJBS1st.txt > " + eachDirAddress + "085log.txt"
        #NJBSline1 = "orthoscopeScripts/Rscript orthoscopeScripts/085_NJBSa.R " + eachDirAddress + "082_trimedBlockExc3rdPhy.txt " + outgroup1 + " K80 " + eachDirAddress + "085_NJBS1st.txt > " + eachDirAddress + "085log.txt"
    elif dataset == "Include3rd":
        #NJBSline1 = "orthoscopeScripts/Rscript orthoscopeScripts/085_NJBSa.R " + eachDirAddress + "082_trimedBlockInc3rdPhy.txt " + outgroup1 + " " + eachDirAddress + "085_NJBS1st.txt > " + eachDirAddress + "085log.txt"
        NJBSline1 = "orthoscopeScripts/Rscript orthoscopeScripts/085_NJBSa.R " + eachDirAddress + "080_trimedCDNAPhy.txt "        + outgroup1 + " TN93 " + eachDirAddress + "085_NJBS1st.txt > " + eachDirAddress + "085log.txt"
    else:
        NJBSline1 = "orthoscopeScripts/fastme -i " + eachDirAddress + "082_trimedAAFastmePhy.txt --protein=WAG -m NJ -b 100 -v 3 -o " + eachDirAddress + "085_NJBS1st.txt > " + eachDirAddress + "085_fastmelog.txt"
else:
    NJBSline1 = "orthoscopeScripts/fastme -i " + eachDirAddress + "082_trimedAAFastmePhy.txt --protein=WAG -m NJ -b 100 -v 3 -o " + eachDirAddress + "085_NJBS1st.txt > " + eachDirAddress + "085_fastmelog.txt"

#print("SSSS22")
#print ("NJBSline1: ", NJBSline1)
subprocess.call(NJBSline1, shell=True)
#exit()



if not os.path.isfile(eachDirAddress + "085_NJBS1st.txt"):
    print("Error1 in gene tree estimation: Cannot estimate the NJ tree.<br>")
    print("The distance matrix was not made.")
    print("Just push Execute again.<br>")
    print("Otherwise, please try followings:<br>")
    print("(1) Exclude species with bad assembly quality.<br>")
    print('(2) Exclude short sequences using "Aligned site rate threshold within unambiguously aligned sites."<br>')
    print("Please note, sometimes bootstrap analysis produces an unavailab alignment for a distance calculation.<br>")
    exit()

NJtreeFile = open(eachDirAddress + "085_NJBS1st.txt")
NJtreeFileCont = list(NJtreeFile)
if not NJtreeFileCont:
    print("Error in gene tree estimation: Cannot estimate the NJ tree.")
    print("Just push Execute again.<br>")
    print("Otherwise, try the DNA analysis if the amino acid analysis got stuck.<br>")
    print("For tree estimation, DNA sequences are better than amino acid sequences in general.<br>")
    exit()

if reconciliation == "treeSearchRearrange":
    #print("### 210_Notung-2.6.jar: 085_NJBS1st.txt.rearrange.0")
    NOTUNG1stLine = "java -jar orthoscopeScripts/210_Notung-2.9.jar -s " + eachDirAddress + "000_uploadedTree.txt -g " + eachDirAddress + "085_NJBS1st.txt --outputdir " + eachDirAddress + " --rearrange --threshold " + str(RearrangementBSthreshold) + " --speciestag prefix  --maxtrees 5 --nolosses --treeoutput nhx > " + eachDirAddress + "085_NOTUNGlog.txt"
    #print("NOTUNG1stLine:", NOTUNG1stLine)
    subprocess.call(NOTUNG1stLine, shell=True)

if reconciliation == "treeSearchRearrange":
    if not os.path.isfile(eachDirAddress + "085_NJBS1st.txt.rearrange.0"):
        print ("Error in NOTUNG analysis: Cannot compare the estimated NJ tree and your species tree.<br>")
        exit()



############## Printout
##print("### makeSummary:100_1stAnalysisSummary.txt")



makeSummary()


treePlotR = "orthoscopeScripts/Rscript orthoscopeScripts/115_treePlotB.R " + eachDirAddress + "100_1stAnalysisSummary.txt " + analysisGroupName + " " + eachDirAddress + "115_1st > " + eachDirAddress + "115_logTreePlotB.txt"
#print("treePlotR: ", treePlotR)
subprocess.call(treePlotR, shell=True)


#print("### cDNAaligHtmlMaker: 230_CDNA.html<br>")
#print("### reorderFas2PhyByTree: 240")
tree4leafOrder = ""
if reconciliation == "treeSearchRearrange":
    tree4leafOrder = "115_1stRearranged_geneTree.txt"
else:
    tree4leafOrder = "085_NJBS1st.txt"


if querySeqChr == "nucl":
    cDNAaligHtmlMaker("054_p2nOutcDNAfas.txt", "052_AA.fas.trm.html", tree4leafOrder, outFile = "230_CDNA.html")


## Raw seqs
reorderDeleteGap_Fas2FasByTree("050_retAAseqs.maf",             tree4leafOrder, outPhyFileName="240_reorderedAAfas.txt")
if querySeqChr == "nucl":
    reorderDeleteGap_Fas2FasByTree("054_p2nOutcDNAfas.txt",         tree4leafOrder, outPhyFileName="240_reorderedCDNAfas.txt")

#reorderFas2PhyByTree("052_AA.fas.trm",               tree4leafOrder, outPhyFileName="240_reorderedTrimedAAPhy.txt")
#reorderFas2PhyByTree("054_p2nOutcDNAfas.txt",        tree4leafOrder, outPhyFileName="240_reorderedCDNAphy.txt")
#rorderPhyFileByTree("080_trimedCDNAPhy.txt",        tree4leafOrder, outPhyFileName="240_reorderedTrimedCDNAphy.txt")
#reorderPhyFileByTree("082_trimedBlockExc3rdPhy.txt", tree4leafOrder, outPhyFileName="240_reorderedTrimedBlockExc3rdPhy.txt")
#reorderPhyFileByTree("082_trimedBlockInc3rdPhy.txt", tree4leafOrder, outPhyFileName="240_reorderedTrimedBlockInc3rdPhy.txt")

#print("### compression: result....zip")
#print("test finished.")
#exit()


change_txt2html()
compression()


#print("### resHtmlMaker: 300_resultsREA.html<br>")
resHtmlMaker()
##print("<br><br>")

### https://qiita.com/fantm21/items/3dc7fbf4e935311488bc
elapsed_time = round((time.time() - startTime),1)

if nums_querySeqs == 1 and int(blastTopHits) == 3:
    print("You used only 1 query sequence.<br>")
    print("Use of 2 or more distantly-related queries is recommended to count orthogroup members.<br>")
    print("Otherwise, select 5 or more for the 'Number of BLAST hits to report per genome' option(below).<br><br>")
print ("Analysis time: {0}".format(elapsed_time) + " seconds")


#if "inouejun-no-Mac-Pro" in EMAIL_SUBJECT_PREFX:
#    dbAddress  = "/dbb/RawDataEns91/"                   ## office
#elif "inouejun-no-MacBook-Pro" in EMAIL_SUBJECT_PREFX:
#    dbAddress   = "/Users/junINOUEpro/db/RawDataEns91/" ## laptop
#else:
#    dbAddress  = "/data/dbfe/RawDataEns91/"             ## fish-evol


#htmlAddress = ""
#if "fish-evol" in EMAIL_SUBJECT_PREFX:
#    htmlAddress = '<br>Finished: <a href="http://fish-evol.unit.oist.jp/orthoscopeWork/' + str(dirname_rand) + '/300_resultsREA.html" target="_blank">result' + str(dirname_rand) + '</a><br><br>'
#elif "sakura" in EMAIL_SUBJECT_PREFX:
#    htmlAddress = '<br>Finished: <a href="http://153.126.167.45/orthoscopeWork/' + str(dirname_rand) + '/300_resultsREA.html" target="_blank">result' + str(dirname_rand) + '</a><br><br>'
#else:
#    htmlAddress = '<br>Finished: <a href="http://localhost/orthoscopeWork/'              + str(dirname_rand) + '/300_resultsREA.html" target="_blank">result' + str(dirname_rand) + '</a><br><br>'
#htmlAddress = make_htmlAddress()
htmlAddress = '<br>Finished: <a href="../orthoscopeWork/' + str(dirname_rand) + '/300_resultsREA.html" target="_blank">result' + str(dirname_rand) + '</a><br><br>'

print(htmlAddress)
#print ("Analysis finished.<br>")


exit()
