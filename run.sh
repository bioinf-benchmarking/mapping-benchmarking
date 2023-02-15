# Runs full pipeline on a clean system with only Snakemake and Conda installed

git checkout master
git pull origin master

pip install -r python_requirements.txt

snakemake --cores 32 --use-conda report.md

# commit report to report-branch, push



