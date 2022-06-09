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

     
