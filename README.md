# hse22_project
project Z-DNA prediction

# Google colab
https://colab.research.google.com/drive/11b1BL_f25FrO5hvJILbK9jBONOEhGlox?usp=sharing 

# Описание выбранных геномов

|*Taxon*|Ascomycetes|||||
|:---|:---|---|---|---|---|
|***Genus***|**Saccharomyces**| **Naumovozyma**||||
|||||||
|***Species***|***Level***|***Size(Mb)***|***GC%***|***Scaffolds***|***Assembly***|
|Saccharomyces arboricola H-6| Chromosome|11,6195|38,851|35|GCA_000292725.1|
|Saccharomyces eubayanus|Chromosome|11,7342|39,9557|24|GCA_001298625.1|
|Saccharomyces paradoxus|Chromosome|12,0928|38,562|17|GCA_002079055.1|
|Naumovozyma castellii CBS 4309| Chromosome|11,2195|36,7559|10|GCA_000237345.1|
|Naumovozyma dairenensis CBS 421|  Chromosome|13,5276|34,1768|11|GCA_000227115.2|  

Все выбранные геномы имеют уровень сборки Chromosome (как рекомендовано в задании), не очень большой размер, а также GC% около 30-40%.
Все данные находятся на Google Drive по [ссылке](https://drive.google.com/drive/folders/1fBILsHN0vko58oxx-LklgGUbbFHdL_J_?usp=sharing).    

# Анализ аннотированных генов  

|***Species***|***Length***|***Annotated genes***|***Annotated genes %***|***Exons %***|
|---|---|---|---|---|
|Saccharomyces arboricola H-6| 11558863 | 5477786 | 0.47390353186122197 | 0.47390353186122197|
|Saccharomyces eubayanus| 11703647 | 8107585 | 0.6927400493196694 | 0.6896107683357162|
|Saccharomyces paradoxus| 12092683 | 8443702 | 0.6982488501517818 | 0.6940445722425701|
|Naumovozyma castellii CBS 4309| 11219539 | 8391917 | 0.7479734238634939 | 0.7439753986326889|
|Naumovozyma dairenensis CBS 421| 13527580 | 8691117 | 0.6424738940741803 | 0.6390060897810251|

### Получение координат аннотированных генов
``` python
def get_gene_coords(file_name):
  genes = {} #словарь
  file1 = open(file_name, "r")
  sum1 = 0
  while True:
      line = file1.readline()
      if not line:
          break
      if ('RefSeq	gene' in line) or ('Genbank	gene' in line):
        stop = int(line.strip().split('\t')[4])
        start = int(line.strip().split('\t')[3])
        gene_id = ((line.strip().split('gene_id')[1]).split(';')[0]).replace('\"', '').replace(' ', '')
        genes[gene_id] = [start, stop]
      
  file1.close
  return genes
```
  
  Таким образом, с помощью этой функции из файла формата .gtf можно получить словарь с ключом - dene_id и координатами начала и конца гена.  
  Пример:  
  ```
  {'NCAS_0A00100': [351, 1271],  
 'NCAS_0A00110': [2003, 2272],  
 'NCAS_0A00120': [3020, 3643],  
 'NCAS_0A00130': [6067, 6759],  
 'NCAS_0A00140': [7200, 10586],  
 ...}
  ```
