# README
先進ゲノム支援 2018年度情報解析中級者講習会の1-1 RNA-seqの各種ツールによる解析のpython実装  
https://github.com/genome-sci/python_bioinfo_2018/tree/master/1-1  
書籍「独習　Pythonバイオ情報解析」の7章「RNA-Seqカウントデータの処理」は上記講習の処理がなされていることが前提となる. 

# 環境
- packages
    ```bash
    sudo apt install g++ \
                     make \
                     unzip \
                     wget \
                     default-jdk \
                     libz-dev  \
                     libpthread-stubs0-dev \
                     libncurses5-dev \
                     ncurses-devel \
                     libbz2-dev \
                     liblzma-dev
    ```


- python version : 3.6
    ```bash
    conda create -n env_name python=3.6
    ```
- python library : bcbio-gff
    ```bash
    pip install bcbio-gff
    ```

# 実行
```bash
python rnaseq.py setup
sudo chmod 777 tools/FastQC/fastqc
python rnaseq.py run
```
