export R_LIBS_USER=/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/wangfang/02.software/miniconda3/envs/R_4.1/lib/R/library

for HAP in 48HAP;do
lambda=1
for c in 200 ; do
for a in 2500;do

ls -d $PWD/*/04.Matrix/FilterMatrix/| grep "${HAP}.*web">count_mtx.filter.list
ls | grep "${HAP}.*web">sample.filter.list
mkdir harmony_${HAP}_c${c}_lambda${lambda}_a${a}_SCT
cd harmony_${HAP}_c${c}_lambda${lambda}_a${a}_SCT
cp ../count_mtx.filter.list ../sample.filter.list  ./

#/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/wangfang/02.software/miniconda3/envs/R_4.1/bin/Rscript 
/hwfssz1/ST_EARTH/P20Z10200N0035/USER/lihuanjin/miniconda3/envs/r-base4/bin/Rscript ../harmony_way_V2_SCT6.R \
-i count_mtx.filter.list \
-s sample.filter.list \
-n 200 \
-a ${a} \
-c ${c} \
-b 20 \
-r 1 \
-U 30 \
-t embryo \
-l ${lambda} \
#-p ^ATM \
#-f ^ATC \
#-q 5 \
#-g 5 \

cd ../
done
done
done
#-i 输入矩阵文件
#-s 输入sample名称文件
#-n 质控，nFeature_RNA > ,默认200
#-a 质控，nFeature_RNA < ,默认3000
#-c 质控，nCount_RNA > ,默认1000
#-b 合并后数据FindNeighbors设置主成分，默认30
#-r 合并后数据设置细胞分群分辨率，默认1
#-U 单个数据RunUMAP,默认30
#-p 计算线粒体或叶绿体基因所占比例,没有可以不设
#-f 计算线粒体或叶绿体基因所占比例,没有可以不设
#-q 质控，删除比例之外的细胞,默认5
#-g 质控，删除比例之外的细胞,默认5

qstat -j $JOB_ID |grep cpu
