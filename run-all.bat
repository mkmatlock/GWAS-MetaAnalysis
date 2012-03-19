python src/metaAnalysis.py --inc-map --verify --update --exclude-traits data\gwas\gwas_trait_exclude.txt --rm-pseudo --output-dir "html_all" 

python src/metaAnalysis.py --inc-map --verify --update --exclude-traits data\gwas\gwas_trait_exclude.txt --rm-pseudo --output-dir "html_partial" --freq 2 

python src/metaAnalysis.py --inc-map --verify --update --exclude-traits data\gwas\gwas_trait_exclude.txt --rm-pseudo --output-dir "html_bustamente" --skip-lists --rm-study "2,3,4,5"
python src/metaAnalysis.py --inc-map --verify --update --exclude-traits data\gwas\gwas_trait_exclude.txt --rm-pseudo --output-dir "html_vamathevan" --skip-lists --rm-study "1,3,4,5"
python src/metaAnalysis.py --inc-map --verify --update --exclude-traits data\gwas\gwas_trait_exclude.txt --rm-pseudo --output-dir "html_kosiol" --skip-lists --rm-study "1,2,4,5"
python src/metaAnalysis.py --inc-map --verify --update --exclude-traits data\gwas\gwas_trait_exclude.txt --rm-pseudo --output-dir "html_bakewell" --skip-lists --rm-study "1,2,3,5"
python src/metaAnalysis.py --inc-map --verify --update --exclude-traits data\gwas\gwas_trait_exclude.txt --rm-pseudo --output-dir "html_nielsen" --skip-lists --rm-study "1,2,3,4"
