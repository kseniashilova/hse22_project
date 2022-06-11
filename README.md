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
# Предсказание участков Z-DNA
Количество предсказанных zhunt участков, их суммарная длина и средняя длина предсказанного участка Z-DNA:  
|***Species***|***Amount***|***Overal length***|***Mean length***| ***Mean ZH-Score*** |
|---|---|---|---|---|
|Saccharomyces arboricola H-6| 7437 | 76850 | 10.33 | 2220.6 |
|Saccharomyces eubayanus| 10197 | 105374 | 10.33 | 2882.8 |
|Saccharomyces paradoxus| 6920 | 69898 | 10.10 | 1769.3 |
|Naumovozyma castellii CBS 4309| 3231 | 32900 | 10.18 | 1491.8 |
|Naumovozyma dairenensis CBS 421| 4617| 49046 | 10.62| 1558.2 |

Файлы с предсказанием Z-DNA сохранены как текстовые файлы в [папке на гугл диске](https://drive.google.com/drive/folders/1fBILsHN0vko58oxx-LklgGUbbFHdL_J_?usp=sharing)  
Гистограммы распределение ZH-Score для пяти геномов:  
  
  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/z1.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/z2.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/z3.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/z4.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/z5.png)
  
    
Как видно из гистограмм, ZH-Score распределен со смещением в сторону небольших значений с некоторыми всплесками (небольшим количеством очень больших значений), но по среднему ZH-Score наиболее подверженен образованию Z-DNA является Saccharomyces eubayanus геном.  
  
  
Затем, создаем .bed файлы с ZH-Score и объединяем предсказанные участки с помощью bedtools merge.  
Все промежуточные файлы для дальнейшней работы находятся в [папке на гугл диске](https://drive.google.com/drive/folders/1fBILsHN0vko58oxx-LklgGUbbFHdL_J_?usp=sharing).  

# Ассоциация предсказанных участков с генами
С помощью ```bedtools slop``` создаем новые .bed файлы, где от координаты гена отступаем по 300bp вправо и влево.  
С помощью ```bedtools intersect``` создаем новые INTERSECTION_\*.bed файлы.  
Для примера я визуализировала первые несколько генов, которые пересеклись с предсказаными Z-DNA для двух геномов *Naumovozyma castellii CBS 4309* и *Saccharomyces eubayanus*.  
Картинка показывает, что Z-DNA действительно накладываются на гены. Значит, предсказания z-hunt и пересечение сработали корректно.  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/vis_gene_example.png)   
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/vis_gene_example2.png)   
Приведу гистограммы, которые демонстрируют долю гена (с учетом 300bp до начала гена), где начинается предсказанная Z-DNA. То есть, на гистограммах большинство значений сосредоточено около нуля, а это означает, что большинство пересечений генов и предсказанных Z-DNA находятся в районе начала генов (промотера), так как доля гена, которая идет ДО Z-DNA близка к нулю.  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/hist_share_Naumovozyma_castellii_CBS_4309.png) ![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/hist_share_Naumovozyma_dairenensis_CBS_421.png) ![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/hist_share_Saccharomyces_arboricola_H-6.png) ![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/hist_share_Saccharomyces_eubayanus.png) ![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/hist_share_Saccharomyces_paradoxus.png)      
Интересно заметить, что на гистограммах есть провалы в районе долей от длины гена, равных 0.25 и 0.75. Эти участки могут соответствовать регионам гена TSS и TES.  
Можно предположить, что Z-DNA мешает транскрипции на этих участках.  

# Определение гомологов. Кластеризация
С помощью ```proteinortho5``` находим гомологичные гены для всех выбранных геномов (предварительно скачиваем .faa файлы с белками). Отбираем среди получившихся кластеров гомолгичных белков те, в которых участвуют все пять геномов.  
Как видно из гистограммы количества геномов в кластере гомологов, кластеров с участием всех пяти геномов наибольшее количество.  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/hists_homology_species.png)    
Затем, нам нужно определить, в каких из кластеров гомологичные гены ассоциированы с предсказанными участками Z-DNA.  
Делаем это с помощью несложного алгоритма (код можно посмотреть в colab, а здесь приведу схему этого алгоритма):  
```
indices = []  
Проходимся по строкам в таблице proteins:  
 flag = True  
 Проходимся по белкам в строке из proteins:  
  
	есть название белка      XP_003673230.1
	ищем название гена в gtf       NCAS_0A02810
	ищем координаты начала гена в словаре     535902	537742
	ищем строку в INTERSECTION _start       (_535902)
		Если такой строки нет, то значит нет Z-DNA
			flag = False - не запоминаем этот кластер
			Пропускаем эту строку в таблице proteins (с гомологичными белками) - BREAK
		Если такая строка есть, 
			Проверяем остальные белки из этой строки (если там список через запятую, то хотя бы один из них) - CONTINUE
 
 if flag:
	indices.append(номер current строки)
```    
```indices``` - это список строк в файле с кластерами (результат работы ```proteinortho5```), где гомологичные гены ассоциированы с Z-DNA из предсказаний z-hunt.    
В результате, если не отдавать предпочтение участкам промотеров в генах, а анализировать гены по всей длине, то подходящих кластеров (где гены ассоциированы с Z-DNA) получилось **1314**.  
Затем, так как нас больше интересуют участки промотеров, сделаем новые файлы slop (только в райноре промотера, отступаем от начала гена 300 влево и 300 вправо, учитывая знак цепи '+' или '-'. Сделаем также новые файлы с помощью ```bedtools intersect```, где будет информация о том, какие из предсказанных Z-DNA находятся в райноре промотера.  
С помощтю алгоритма, схема которого представлена выше, найдем новые кластеры, где гомологичные гены пересекаются с предсказаниями z-hunt в районе *промотера*.  
Таких кластеров получилось всего лишь **42**.  
Запишем среднее значение Z-Hunt-Score для этих кластеров и возьмем первые 10 кластеров с наибольшим средним ZH-Score. Получаем следующую таблицу (*str_number* - номер кластера в файле):  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/table_zhscores.PNG)  
Как мы видим, наибольшим средний ZH Score сильно отличается от других (он около 35000). Это означает что на первый из этого списка кластер нужно обратить особое внимание при анализе.   
  
  Затем получаем таблицу, где запоминаем координаты генов, которые соответствуют белкам и выбранных десяти кластеров, а также запоминаем координаты участков Z-DNA, которые пересекаются с промотерами генов.  
  |index|\#Species|Genes|Alg\.-Conn.|GCA\_000227115.2.faa|GCA\_000237345.1.faa|GCA\_000292725.1.faa|GCA\_001298625.1.faa|GCA\_002079055.1.faa|ZH-Scores|Mean-ZH-Score|Gene\_GCA_000237345.1.faa|Gene\_GCA_000227115.2.faa|Gene\_GCA_000292725.1.faa|Gene\_GCA_001298625.1.faa|Gene\_GCA_002079055.1.faa|ZDNA\_GCA_000237345.1.faa|ZDNA\_GCA_000227115.2.faa|ZDNA\_GCA_000292725.1.faa|ZDNA\_GCA_001298625.1.faa|ZDNA\_GCA_002079055.1.faa|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|0|5|5|1\.0|XP\_003672405.1|XP\_003678511.1|EJS42635\.1|XP\_018220820.1|XP\_033768165.1|\[3428.529,883.5764,2752.447,138924.1,28780.5]|34953\.830480000004|371235,372539|676541,677836|944561,945835|975968,977242|961470,962804|371079,371091|676334,676340|944275,944287|975754,975775|961285,961303|
|1|5|5|1\.0|XP\_003671415.2|XP\_003675812.1|EJS41822\.1|XP\_018221672.1|XP\_033769035.1|\[780.9275,766.6232,612.3848,783.823,28780.5]|6344\.8517|936565,937947|941809,943185|54008,55432|58306,59742|74401,75831|936291,936303|941509,941517|53814,53824|58420,58430|74237,74256|
|2|5|5|1\.0|XP\_003668882.1|XP\_003677216.1|EJS44048\.1|XP\_018222699.1|XP\_033765848.1|\[826.2209,883.5764,3855.439,904.32,23787.35]|6051\.38126|763919,764569|1475362,1476009|199998,200588|224756,225310|218622,219179|763681,763692|1475154,1475162|200274,200287|224518,224528|218611,218628|
|3|5|5|1\.0|XP\_003670308.1|XP\_003675532.1|EJS41431\.1|XP\_018218919.1|XP\_033769652.1|\[1280.855,16271.57,900.2643,883.5764,2924.216]|4452\.09634|324446,326548|512966,515020|311450,313543|310011,312110|314358,316451|324373,324385|512831,512856|311199,311211|310221,310229|314150,314171|
|4|5|5|1\.0|XP\_003671627.1|XP\_003673710.1|EJS44026\.1|XP\_018222664.1|XP\_033765814.1|\[905.6763,16271.57,785.5658,752.603,766.6232]|3896\.40766|1533902,1534597|512715,513416|134166,134861|158574,159263|154568,155257|1533630,1533642|512831,512856|133936,133946|158274,158278|154601,154609|
|5|5|5|1\.0|XP\_003672661.1|XP\_003673485.1|EJS41221\.1|XP\_018219904.1|XP\_033768967.1|\[1030.205,16271.57,650.9198,908.3955,583.4285]|3888\.90376|1075589,1075918|512560,512892|911,1207|687217,687513|683517,683813|1075759,1075773|512831,512856|91428,91436|687390,687400|683220,683232|
|6|5|5|1\.0|XP\_003668512.1|XP\_003674949.1|EJS43948\.1|XP\_018222841.1|XP\_033765995.1|\[915.9191,883.5764,1820.585,2091.083,12723.6]|3686\.9527|905447,906721|587643,589010|510373,511635|545348,546610|537305,538567|905332,905344|587897,587903|510104,510121|545066,545079|537024,537046|
|7|5|5|1\.0|XP\_003671326.2|XP\_003675727.1|EJS43925\.1|XP\_018222504.1|XP\_033766031.1|\[752.603,1820.585,1232.502,13713.99,766.6232]|3657\.26064|748210,752190|712935,716999|65868,69923|39614,43672|79183,83238|747959,747973|712707,712726|65598,65615|39432,39447|79428,79438|
|8|5|5|1\.0|XP\_003671423.2|XP\_003675821.1|EJS41828\.1|XP\_018221680.1|XP\_033769044.1|\[3488.771,783.823,10894.72,883.5764,790.8198]|3368\.3420399999995|956018,957778|959388,961136|71823,73544|75824,77593|92065,93822|955823,955838|959281,959291|72044,72058|75659,75665|91978,91992|
|9|5|5|1\.0|XP\_003670078.1|XP\_003674480.1|EJS42479\.1|XP\_018219947.1|XP\_033768194.1|\[980.8116,705.4245,650.9198,4860.206,8922.646]|3224\.00158|15368,17134|15301,17109|11481,13226|27240,28985|21970,23715|15546,15558|15030,15042|11369,11377|27077,27091|21762,21788|    
  
    
      
# Визуализация десяти кластеров *DNA Features Viewer*  (ZH-Scores обозначены в заголовке картинок).   
В [папке](https://github.com/kseniashilova/hse22_project/tree/main/clusts) находятся картинки с визуализацией пересечения Z-DNA и генов из кластера, а также файлы в формате .faa для множественного выравнивания.  
### Визуализация кластеров и ассоциированных предсказанных Z-DNA с пересечением с генами в районе промотера:  
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust0.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust1.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust2.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust3.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust4.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust5.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust6.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust7.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust8.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust9.png)


  
  Хотя на некоторых картинках и кажется, что в гомологичных генах координаты одинаковы, но на самом деле они отличаются. Это можно увидеть в таблице выше (10 столбцов в конце содержат координаты генов и Z-DNA, по которым строилась визуализация пересечения).    
    
   
Рассмотрим похожую визуализацию, но со всеми предсказанными Z-DNA (Не только на пересечении с промотером):   
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust0_other_ZDNA.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust1_other_ZDNA.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust2_other_ZDNA.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust3_other_ZDNA.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust4_other_ZDNA.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust5_other_ZDNA.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust6_other_ZDNA.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust7_other_ZDNA.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust8_other_ZDNA.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/clusts/clust9_other_ZDNA.png)

     
# Множественное выравнивание с помощью ClustalW (Mega-X)  
С помощью ClustlW алгоритма (на десктопной MEGA-X) выравниваем последовательности кластеров. В [папке](https://github.com/kseniashilova/hse22_project/tree/main/alignment) можно найти файлы с выравниваниями.  
Получились приблизительно такие картинки (для контроля того, что выравнивание получилось корректным) - Начало выравненных последовательностей и середина.    
**Clust 0**  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/0_start.PNG)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/0_mid.PNG)    
**Clust 1**  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/1_start.PNG)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/1_mid.PNG)  
**Clust 2**  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/2_start.PNG)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/2_mid.PNG)  
**Clust 3**  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/3_start.PNG)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/3_mid.PNG)    
**Clust 4**  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/4_start.PNG)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/4_mid.PNG)    
**Clust 5**  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/5_start.PNG)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/5_mid.PNG)    
**Clust 6**  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/6_start.PNG)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/6_mid.PNG)    
**Clust 7**  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/7_start.PNG)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/7_mid.PNG)    
**Clust 8**  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/8_start.PNG)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/8_mid.PNG)    
**Clust 9**  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/9_start.PNG)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/alignment_pictures/9_mid.PNG)
  
    
Как можно видеть на картинках выше, выравнивание прошло успешно. В принципе, несмотря на то, что различий в аминокислотах много, похожие участки все равно хорошо прослеживаются.  

  
# Функции генов из кластеров  
Таблица ниже показывает ID генов из кластеров (5 правых столбцов)  
|index|\#Species|Genes|Alg\.-Conn.|GCA\_000227115.2.faa|GCA\_000237345.1.faa|GCA\_000292725.1.faa|GCA\_001298625.1.faa|GCA\_002079055.1.faa|ZH-Scores|Mean-ZH-Score|GeneID\_GCA_000237345.1.faa|GeneID\_GCA_000227115.2.faa|GeneID\_GCA_000292725.1.faa|GeneID\_GCA_001298625.1.faa|GeneID\_GCA_002079055.1.faa|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|0|5|5|1\.0|XP\_003672405.1|XP\_003678511.1|EJS42635\.1|XP\_018220820.1|XP\_033768165.1|\[3428.529,883.5764,2752.447,138924.1,28780.5]|34953\.830480000004|NCAS\_0J01940|NDAI\_0J02700|SU7\_2398|DI49\_3933|SPAR\_L04410|
|1|5|5|1\.0|XP\_003671415.2|XP\_003675812.1|EJS41822\.1|XP\_018221672.1|XP\_033769035.1|\[780.9275,766.6232,612.3848,783.823,28780.5]|6344\.8517|NCAS\_0C04580|NDAI\_0G03950|SU7\_3101|DI49\_2282|SPAR\_O00310|
|2|5|5|1\.0|XP\_003668882.1|XP\_003677216.1|EJS44048\.1|XP\_018222699.1|XP\_033765848.1|\[826.2209,883.5764,3855.439,904.32,23787.35]|6051\.38126|NCAS\_0F03790|NDAI\_0B06070|SU7\_0859|DI49\_1433|SPAR\_E01050|
|3|5|5|1\.0|XP\_003670308.1|XP\_003675532.1|EJS41431\.1|XP\_018218919.1|XP\_033769652.1|\[1280.855,16271.57,900.2643,883.5764,2924.216]|4452\.09634|NCAS\_0C01760|NDAI\_0E02480|SU7\_3544|DI49\_5347|SPAR\_P01550|
|4|5|5|1\.0|XP\_003671627.1|XP\_003673710.1|EJS44026\.1|XP\_018222664.1|XP\_033765814.1|\[905.6763,16271.57,785.5658,752.603,766.6232]|3896\.40766|NCAS\_0A07710|NDAI\_0H02100|SU7\_0837|DI49\_1396|SPAR\_E00710|
|5|5|5|1\.0|XP\_003672661.1|XP\_003673485.1|EJS41221\.1|XP\_018219904.1|XP\_033768967.1|\[1030.205,16271.57,650.9198,908.3955,583.4285]|3888\.90376|NCAS\_0A05440|NDAI\_0K02270|SU7\_u0006|DI49\_4768|SPAR\_N03370|
|6|5|5|1\.0|XP\_003668512.1|XP\_003674949.1|EJS43948\.1|XP\_018222841.1|XP\_033765995.1|\[915.9191,883.5764,1820.585,2091.083,12723.6]|3686\.9527|NCAS\_0B04930|NDAI\_0B02340|SU7\_0963|DI49\_1587|SPAR\_E02520|
|7|5|5|1\.0|XP\_003671326.2|XP\_003675727.1|EJS43925\.1|XP\_018222504.1|XP\_033766031.1|\[752.603,1820.585,1232.502,13713.99,766.6232]|3657\.26064|NCAS\_0C03720|NDAI\_0G03060|SU7\_0976|DI49\_1620|SPAR\_F00270|
|8|5|5|1\.0|XP\_003671423.2|XP\_003675821.1|EJS41828\.1|XP\_018221680.1|XP\_033769044.1|\[3488.771,783.823,10894.72,883.5764,790.8198]|3368\.3420399999995|NCAS\_0C04670|NDAI\_0G04030|SU7\_3107|DI49\_2290|SPAR\_O00400|
|9|5|5|1\.0|XP\_003670078.1|XP\_003674480.1|EJS42479\.1|XP\_018219947.1|XP\_033768194.1|\[980.8116,705.4245,650.9198,4860.206,8922.646]|3224\.00158|NCAS\_0B00190|NDAI\_0E00190|SU7\_2413|DI49\_3966|SPAR\_M00070|      
      
          
	  
Добавим функции белков в таблицу (столбцы справа) из файлов .gbff:    
     
     
|index|\#Species|Genes|Alg\.-Conn.|GCA\_000227115.2.faa|GCA\_000237345.1.faa|GCA\_000292725.1.faa|GCA\_001298625.1.faa|GCA\_002079055.1.faa|ZH-Scores|Mean-ZH-Score|GeneID\_GCA_000237345.1.faa|GeneID\_GCA_000227115.2.faa|GeneID\_GCA_000292725.1.faa|GeneID\_GCA_001298625.1.faa|GeneID\_GCA_002079055.1.faa|Note\_GCA_000237345.1.faa|Note\_GCA_000227115.2.faa|Note\_GCA_000292725.1.faa|Note\_GCA_001298625.1.faa|Note\_GCA_002079055.1.faa|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|0|5|5|1\.0|XP\_003672405.1|XP\_003678511.1|EJS42635\.1|XP\_018220820.1|XP\_033768165.1|\[3428.529,883.5764,2752.447,138924.1,28780.5]|34953\.830480000004|NCAS\_0J01940|NDAI\_0J02700|SU7\_2398|DI49\_3933|SPAR\_L04410|                     ancestral locus Anc\_4.317|                     similar to Saccharomyces cerevisiae CAR2 \(YLR438W);                     ancestral locus Anc_4.317|                     similar to Saccharomyces cerevisiae CAR2 \(YLR438W)|                     ancestor homolog: Anc\_4.317; saccharomyces                     cerevisiae ortholog: YLR438W|                     L-ornithine transaminase \(OTAse);                     similar to YLR438W|
|1|5|5|1\.0|XP\_003671415.2|XP\_003675812.1|EJS41822\.1|XP\_018221672.1|XP\_033769035.1|\[780.9275,766.6232,612.3848,783.823,28780.5]|6344\.8517|NCAS\_0C04580|NDAI\_0G03950|SU7\_3101|DI49\_2282|SPAR\_O00310|                     ancestral locus Anc\_3.46|                     similar to Saccharomyces cerevisiae TRM13                     \(YOL125W); ancestral locus Anc_3.46|                     similar to Saccharomyces cerevisiae YOL125W|                     ancestor homolog: Anc\_3.46; saccharomyces                     cerevisiae ortholog: YOL125W|                     2'-O-methyltransferase;                     similar to YOL125W|
|2|5|5|1\.0|XP\_003668882.1|XP\_003677216.1|EJS44048\.1|XP\_018222699.1|XP\_033765848.1|\[826.2209,883.5764,3855.439,904.32,23787.35]|6051\.38126|NCAS\_0F03790|NDAI\_0B06070|SU7\_0859|DI49\_1433|SPAR\_E01050|                     ancestral locus Anc\_3.532|                     similar to Saccharomyces cerevisiae YER034W;                     ancestral locus Anc\_3.532|                     similar to Saccharomyces cerevisiae YER034W|                     ancestor homolog: Anc\_3.532; saccharomyces                     cerevisiae ortholog: YER034W|                     similar to YER034W|
|3|5|5|1\.0|XP\_003670308.1|XP\_003675532.1|EJS41431\.1|XP\_018218919.1|XP\_033769652.1|\[1280.855,16271.57,900.2643,883.5764,2924.216]|4452\.09634|NCAS\_0C01760|NDAI\_0E02480|SU7\_3544|DI49\_5347|SPAR\_P01550|                     ancestral locus Anc\_8.607|                     similar to Saccharomyces cerevisiae HOS3 \(YPL116W);                     ancestral locus Anc_8.607|                     similar to Saccharomyces cerevisiae HOS3 \(YPL116W)|                     ancestor homolog: Anc\_8.607; saccharomyces                     cerevisiae ortholog: YPL116W|                     Trichostatin A-insensitive homodimeric histone                     deacetylase \(HDAC);                     similar to YPL116W|
|4|5|5|1\.0|XP\_003671627.1|XP\_003673710.1|EJS44026\.1|XP\_018222664.1|XP\_033765814.1|\[905.6763,16271.57,785.5658,752.603,766.6232]|3896\.40766|NCAS\_0A07710|NDAI\_0H02100|SU7\_0837|DI49\_1396|SPAR\_E00710|                     ancestral locus Anc\_7.146|                     similar to Saccharomyces cerevisiae NOP16                     \(YER002W); ancestral locus Anc_7.146|                     similar to Saccharomyces cerevisiae NOP16                     \(YER002W)|                     ancestor homolog: Anc\_7.146; saccharomyces                     cerevisiae ortholog: YER002W|                     Constituent of 66S pre-ribosomal particles;                     similar to YER002W|
|5|5|5|1\.0|XP\_003672661.1|XP\_003673485.1|EJS41221\.1|XP\_018219904.1|XP\_033768967.1|\[1030.205,16271.57,650.9198,908.3955,583.4285]|3888\.90376|NCAS\_0A05440|NDAI\_0K02270|SU7\_u0006|DI49\_4768|SPAR\_N03370|                     ancestral locus Anc\_6.348|                     similar to Saccharomyces cerevisiae YNR034W-A and                     YCR075W-A; ancestral locus Anc\_6.348|                     similar to Saccharomyces cerevisiae YNR034W-A|                     ancestor homolog: Anc\_6.348; saccharomyces                     cerevisiae ortholog: YNR034W-A|                     Protein with a possible role in tRNA export;                     similar to YNR034W|
|6|5|5|1\.0|XP\_003668512.1|XP\_003674949.1|EJS43948\.1|XP\_018222841.1|XP\_033765995.1|\[915.9191,883.5764,1820.585,2091.083,12723.6]|3686\.9527|NCAS\_0B04930|NDAI\_0B02340|SU7\_0963|DI49\_1587|SPAR\_E02520|                     ancestral locus Anc\_8.246|                     similar to Saccharomyces cerevisiae PDA1 \(YER178W);                     ancestral locus Anc_8.246|                     similar to Saccharomyces cerevisiae PDA1 \(YER178W)|                     ancestor homolog: Anc\_8.246; saccharomyces                     cerevisiae ortholog: YER178W|                     E1 alpha subunit of the pyruvate dehydrogenase                     \(PDH) complex;                     similar to YER178W|
|7|5|5|1\.0|XP\_003671326.2|XP\_003675727.1|EJS43925\.1|XP\_018222504.1|XP\_033766031.1|\[752.603,1820.585,1232.502,13713.99,766.6232]|3657\.26064|NCAS\_0C03720|NDAI\_0G03060|SU7\_0976|DI49\_1620|SPAR\_F00270|                     ancestral locus Anc\_8.28|                     similar to Saccharomyces cerevisiae RPO41                     \(YFL036W); ancestral locus Anc_8.28|                     similar to Saccharomyces cerevisiae RPO41                     \(YFL036W)|                     ancestor homolog: Anc\_8.28; saccharomyces                     cerevisiae ortholog: YFL036W|                     Mitochondrial RNA polymerase;                     similar to YFL036W|
|8|5|5|1\.0|XP\_003671423.2|XP\_003675821.1|EJS41828\.1|XP\_018221680.1|XP\_033769044.1|\[3488.771,783.823,10894.72,883.5764,790.8198]|3368\.3420399999995|NCAS\_0C04670|NDAI\_0G04030|SU7\_3107|DI49\_2290|SPAR\_O00400|                     ancestral locus Anc\_3.58|                     similar to Saccharomyces cerevisiae TRF5 \(YNL299W)                     and PAP2 (YOL115W); ancestral locus Anc_3.58|                     similar to Saccharomyces cerevisiae TRF4 \(YOL115W)|                     ancestor homolog: Anc\_3.58; saccharomyces                     cerevisiae ortholog: YOL115W|                     Non-canonical poly\(A) polymerase;                     similar to YOL115W|
|9|5|5|1\.0|XP\_003670078.1|XP\_003674480.1|EJS42479\.1|XP\_018219947.1|XP\_033768194.1|\[980.8116,705.4245,650.9198,4860.206,8922.646]|3224\.00158|NCAS\_0B00190|NDAI\_0E00190|SU7\_2413|DI49\_3966|SPAR\_M00070|                     ancestral locus Anc\_8.866|                     similar to Saccharomyces cerevisiae RSC9 \(YML127W);                     ancestral locus Anc_8.866|                     similar to Saccharomyces cerevisiae RSC9 \(YML127W)|                     ancestor homolog: Anc\_8.866; saccharomyces                     cerevisiae ortholog: YML127W|                     Component of the RSC chromatin remodeling complex;                     similar to YML127W|
    
        
Скриншот функций (для удобства чтения):  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/functions.png)
	  
	  
**Во втором геноме:**  
```YLR438W``` - L-ornithine transaminase (OTAse); catalyzes the second step of arginine degradation, expression is dually-regulated by allophanate induction and a specific arginine induction process; not nitrogen catabolite repression sensitive; protein abundance increases in response to DNA replication stress; human homolog OAT complements yeast null mutan.  
```YOL125W``` - 2'-O-methyltransferase; responsible for modification of tRNA at position 4; C-terminal domain has similarity to Rossmann-fold (RFM) superfamily of RNA methyltransferases.  
```YER034W``` - hypothetical protein; non-essential gene; expression induced upon calcium shortage; protein abundance increases in response to DNA replication stress.  ```YPL116W``` - Trichostatin A-insensitive homodimeric histone deacetylase (HDAC); specificity in vitro for histones H3, H4, H2A, and H2B; similar to Hda1p, Rpd3p, Hos1p, and Hos2p; deletion results in increased histone acetylation at rDNA repeats.  
```YER002W``` - Constituent of 66S pre-ribosomal particles; involved in 60S ribosomal subunit biogenesis.  
```YNR034W-A and YCR075W-A``` - hypothetical protein; expressed during diauxic shift and stationary phase, and negatively regulated by glucose; expression is regulated by Msn2p/Msn4p; overexpression slows down progression through meiosis and improves fermentative efficiency; YNR034W-A has a paralog, YCR075W-A, that arose from the whole genome duplication. Subunit of the EGO/GSE complex; the vacuolar/endosomal membrane associated EGO/GSE complex regulates exit from rapamycin-induced growth arrest, stimulating microautophagy and sorting of Gap1p from the endosome to the plasma membrane; identified by homology to Ashbya gossypii; EGO2 has a paralog, EGO4, that arose from the whole genome duplication.  
```YER178W``` - E1 alpha subunit of the pyruvate dehydrogenase (PDH) complex; catalyzes the direct oxidative decarboxylation of pyruvate to acetyl-CoA; phosphorylated; regulated by glucose; PDH complex is concentrated in spots within the mitochondrial matrix, often near the ERMES complex and near peroxisomes.  
```YFL036W``` - Mitochondrial RNA polymerase; single subunit enzyme similar to those of T3 and T7 bacteriophages; requires a specificity subunit encoded by MTF1 for promoter recognition; Mtf1p interacts with and stabilizes the Rpo41p-promoter complex, enhancing DNA bending and melting to facilitate pre-initiation open complex formation; Rpo41p also synthesizes RNA primers for mitochondrial DNA replication.  
```YOL115W``` - Non-canonical poly(A) polymerase; involved in nuclear RNA degradation as a component of TRAMP; catalyzes polyadenylation of hypomodified tRNAs, and snoRNA and rRNA precursors; required for mRNA surveillance and maintenance of genome integrity, serving as a link between RNA and DNA metabolism; overlapping but non-redundant functions with Trf5p; relocalizes to cytosol in response to hypoxia.  
```YML127W``` - Component of the RSC chromatin remodeling complex; DNA-binding protein involved in the synthesis of rRNA and in transcriptional repression and activation of genes regulated by the Target of Rapamycin (TOR) pathway.   
  
**В третьем геноме:** (Кластер 8) ```YOL115W``` - Non-canonical poly(A) polymerase; involved in nuclear RNA degradation as a component of TRAMP; catalyzes polyadenylation of hypomodified tRNAs, and snoRNA and rRNA precursors; required for mRNA surveillance and maintenance of genome integrity, serving as a link between RNA and DNA metabolism; overlapping but non-redundant functions with Trf5p; relocalizes to cytosol in response to hypoxia.  
  
  
Все остальные замечания про 'similar to smth' одинаковы для гомологичных генов.

  
# G-квадруплексы   
  
С помощью программы ```pqsfinder``` обработаем последовательности fasta файлов для тех же геномов и проанализируем предсказанные G-квадруплексы.  
Код запуска ```pqsfinder``` на **R** можно посмотреть в [гугл ноутбуке](https://colab.research.google.com/drive/1vkTcgQBPDg0ZHQJTi3rgFe-JoBdgyr_Z?usp=sharing). Файлы, полученных с помощью этого кода лежат в этом репозитории в [папке quadruplex](https://github.com/kseniashilova/hse22_project/tree/main/quadruplex).  
Анализ можно посмотреть в одном из последних разделов [ноутбука на Python](https://colab.research.google.com/drive/11b1BL_f25FrO5hvJILbK9jBONOEhGlox?usp=sharing) (все же в Python привычнее и удобнее делать анализ файлов и строить таблицы и графики).   
Таблица ниже показывает количество предсказанных G-квадруплексов, их суммарную длину, среднюю длину и средний Score (из результата работы ```pqsfinder```).  

|***Species***|***Amount of G-quads***|***Overal length***|***Mean length***|***Mean Score***|
|:---|:---|---|---|---|
|Saccharomyces arboricola H-6|74|1597|21.58|67.3|
|Saccharomyces eubayanus|115|2682|23.3|63.2|
|Saccharomyces paradoxus|110|2523|22.9|64.7|
|Naumovozyma castellii CBS 4309| 50 | 1145 | 22.9| 65.32|
|Naumovozyma dairenensis CBS 421| 176 |3528| 20.05|68.6|
  
  
Графики показывают распределение Score из результата предсказания квадруплексов для пяти геномов.  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/hist_score_quad_Naumovozyma%20castellii%20CBS%204309.png) 
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/hist_score_quad_Naumovozyma%20dairenensis%20CBS%20421.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/hist_score_quad_Saccharomyces%20arboricola%20H-6.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/hist_score_quad_Saccharomyces%20eubayanus.png)
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/hist_score_quad_Saccharomyces%20paradoxus.png)    
    
      
 Пример визуализации первых десяти пересечений для двух геномов (*Saccharomyces_eubayanus* и *Naumovozyma_castellii_CBS_4309*):  
 ![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/quad_example.png)  
 ![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/quad_example2.png)    
     
  Гистограммы ниже показывают долю гена, которая идет ДО момента старта предсказанного квадруплекса. То есть, из-за того, что по горизонтальной оси располагается доля гена (координата от старта, нормированная на длину), то и значения там от 0 до 1. По вертикальной оси количество предсказанных G-квадруплексов, которые начинаются на определенной координате от начала гена. На гистограммах заметно, что большинство квадруплексов начинаются в начале гена (доля длины ДО квадруплекса близка к нулю - самый высокий столбец на некоторых диаграммах), а это означает, что и большинство квадруплексов располагаются в районе промотера.    
  ![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/quad_hist_share_Naumovozyma_castellii_CBS_4309.png) 
  ![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/quad_hist_share_Naumovozyma_dairenensis_CBS_421.png) 
  ![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/quad_hist_share_Saccharomyces_arboricola_H-6.png) 
  ![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/quad_hist_share_Saccharomyces_eubayanus.png) ![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/quad_hist_share_Saccharomyces_paradoxus.png)    
  Нужно заметить, что для первого генома (и для третьего в какой-то степени) плозо соблюдается правило всплекса в райноре промотера. Для них наибольшее число предсказанных G-квадруплексов находятся где-то в середине генов (могут быть совершенно разные участки, для более качественного анализа нужно посмотреть на аннотации генов как в одном из предыдущих ДЗ).    
**Консервативность квадруплексов в кластерах гомологичных генов**:  
Аналогично обработке кластеров (с поиком консервативности Z-DNA) обрабатываем файл с предсказаниями для G-квадруплексов.  
Если рассматривать кластеры, где участвуют все 5 геномов, то кластера, где сохраняется консервативность квдруплекса, всего два.      
|index|\# Species|Genes|Alg\.-Conn.|GCA\_000227115.2.faa|GCA\_000237345.1.faa|GCA\_000292725.1.faa|GCA\_001298625.1.faa|GCA\_002079055.1.faa|Scores|
|---|---|---|---|---|---|---|---|---|---|
|1809|5|5|1\.0|XP\_003670050.1|XP\_003673246.1|EJS42813\.1|XP\_018220399.1|XP\_033767742.1|97,73,73,63,63|
|2384|5|5|1\.0|XP\_003670864.1|XP\_003677464.1|EJS42009\.1|XP\_018219741.1|XP\_033768806.1|64,73,64,70,52|  
  
  
При добавлении других кластеров (с участием только 4 из 5 геномов) тоже остаются подходящими только 2 кластера из 5507, которые сформировала программа ```proteinortho5```.  Значит, будем анализировать два этих кластера.  
Визуализация расположения квадруплексов на генах из выбранных кластеров:  
![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/vis_clust0.png) ![](https://github.com/kseniashilova/hse22_project/blob/main/pictures/vis_clust1.png)
**Выводы**
