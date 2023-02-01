import gzip
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from ftplib import FTP_TLS

from BCBio import GFF


class Scrna_seq:
    def __init__(self, fastq_dir, reference_dir, 
                 fastqc_dir_1, fastqc_dir_2, 
                 trimming_dir, hisat_dir, 
                 featurecounts_dir, tool_dir, 
                 sratoolkit_ver="3.0.1", 
                 fastqc_ver="0.11.9", 
                 trimmomatic_ver="0.39", 
                 hisat_ver="2.2.1", 
                 subread_ver="2.0.2", 
                 samtools_ver="1.16.1", 
                 threads=1, 
                 reference_fna="s288c.fna", 
                 reference_gff="s288c.gff"):
        self.fastq_dir = fastq_dir
        self.reference_dir = reference_dir
        self.fastqc_dir_1 = fastqc_dir_1
        self.fastqc_dir_2 = fastqc_dir_2
        self.trimming_dir = trimming_dir
        self.hisat_dir = hisat_dir
        self.featurecounts_dir = featurecounts_dir
        self.tool_dir = tool_dir

        self.sratoolkit_ver = sratoolkit_ver
        self.fastqc_ver = fastqc_ver
        self.trimmomatic_ver = trimmomatic_ver
        self.hisat_ver = hisat_ver
        self.subread_ver = subread_ver
        self.samtools_ver = samtools_ver

        self.fastq_dump = os.path.join(
            self.tool_dir, f"sratoolkit.{self.sratoolkit_ver}-ubuntu64", 
            "bin", "fastq-dump")
        self.trimmomatic = os.path.join(
            self.tool_dir, 
            f"Trimmomatic-{self.trimmomatic_ver}", 
            f"trimmomatic-{self.trimmomatic_ver}.jar")
        self.fastqc = os.path.join(self.tool_dir, "FastQC", "fastqc")
        self.hisat2 = os.path.join(self.tool_dir, 
                                   f"hisat2-{self.hisat_ver}", "hisat2")
        self.hisat2_build = os.path.join(self.tool_dir, 
                                         f"hisat2-{self.hisat_ver}", "hisat2-build")
        self.featureCounts = os.path.join(
            self.tool_dir, 
            f"subread-{self.subread_ver}", 
            "bin", "featureCounts")
        self.samtools = os.path.join(self.tool_dir, 
                                     f"samtools-{self.samtools_ver}", "samtools")
        
        self.reference_fna = os.path.join(self.reference_dir, reference_fna)
        self.reference_gff = os.path.join(self.reference_dir, reference_gff)
        self.threads = str(threads)
        
        os.makedirs(self.fastq_dir, exist_ok=True)
        os.makedirs(self.reference_dir, exist_ok=True)
        os.makedirs(self.fastqc_dir_1, exist_ok=True)
        os.makedirs(self.fastqc_dir_2, exist_ok=True)
        os.makedirs(self.trimming_dir, exist_ok=True)
        os.makedirs(self.hisat_dir, exist_ok=True)
        os.makedirs(self.featurecounts_dir, exist_ok=True)
        os.makedirs(self.tool_dir, exist_ok=True)
        
        trimmomatic_dir = os.path.dirname(self.trimmomatic)
        self.trim_adapt_dir = os.path.join(trimmomatic_dir, "adapters")


    def download_software(self, file_url, save_file, force_download, 
                          make=False, makedir=None, makefile=None, 
                          configure=False):
        save_path = os.path.join(self.tool_dir, save_file)

        if (not os.path.exists(save_path)) or force_download:
            print(f"Downloading {save_path}")
            urllib.request.urlretrieve(file_url, save_path)
            
            if save_file.endswith(".zip"):
                shutil.unpack_archive(save_path, self.tool_dir)
            elif save_file.endswith(".tar.gz"):
                with tarfile.open(save_path, "r:gz") as tar:
                    tar.extractall(path=self.tool_dir)
            elif save_file.endswith(".tar.bz2"):
                with tarfile.open(save_path, "r:bz2") as tar:
                    tar.extractall(path=self.tool_dir)

            if make:
                if configure:
                    subprocess.run(["./configure"], cwd=makedir)

                if makefile is not None:
                    cmd = ["make", "-f", makefile]
                else:
                    cmd = ["make"]
                subprocess.run(cmd, cwd=makedir)


    def setup_software(self, force_download=False):
        # Download sratoolkit
        save_file = f"sratoolkit.{self.sratoolkit_ver}-ubuntu64.tar.gz"
        file_url = f"ftp://ftp.ncbi.nlm.nih.gov/sra/sdk/{self.sratoolkit_ver}/{save_file}"
        self.download_software(file_url, save_file, force_download)

        # # Download fastqc
        save_file = f"fastqc_v{self.fastqc_ver}.zip"
        file_url = f"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/{save_file}"
        self.download_software(file_url, save_file, force_download)

        # Download trimmomatic
        save_file = f"Trimmomatic-{self.trimmomatic_ver}.zip"
        file_url = f"http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/{save_file}"
        self.download_software(file_url, save_file, force_download)

        # Download hisat2
        save_file = f"v{self.hisat_ver}.tar.gz"
        file_url = f"https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/{save_file}"
        makedir = os.path.join(self.tool_dir, f"hisat2-{self.hisat_ver}")
        self.download_software(file_url, save_file, force_download, make=True, makedir=makedir)

        # Download subread
        save_file = f"{self.subread_ver}.tar.gz"
        file_url = f"https://github.com/ShiLab-Bioinformatics/subread/archive/refs/tags/{save_file}"
        makedir = os.path.join(self.tool_dir, f"subread-{self.subread_ver}", "src")
        self.download_software(file_url, save_file, force_download, 
                               make=True, makedir=makedir, makefile="Makefile.Linux")
        
        # Download samtools
        save_file = f"samtools-{self.samtools_ver}.tar.bz2"
        file_url = f"https://github.com/samtools/samtools/releases/download/1.16.1/{save_file}"
        makedir = os.path.join(self.tool_dir, f"samtools-{self.samtools_ver}")
        self.download_software(file_url, save_file, force_download, 
                               make=True, makedir=makedir, configure=True)


    def download_data(self, force_download=False):
        sample_urls = [
            "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX135/SRX135198/SRR453566/SRR453566.sra", 
            "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX135/SRX135198/SRR453567/SRR453567.sra", 
            "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX135/SRX135198/SRR453568/SRR453568.sra", 
            "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX135/SRX135710/SRR453569/SRR453569.sra", 
            "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX135/SRX135710/SRR453570/SRR453570.sra", 
            "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX135/SRX135710/SRR453571/SRR453571.sra"
        ]
        for sample_url in sample_urls:
            filename = os.path.basename(sample_url)
            save_path = os.path.join(self.fastq_dir, filename)

            if (not os.path.exists(save_path)) or force_download:
                print(f"Downloading {save_path}")
                urllib.request.urlretrieve(sample_url, save_path)
                self.dump(save_path)
        
        # ncbiのftpからダウンロードするとなぜか文字化けが発生してしまうので、
        # 仕方なく講義資料からダウンロード
        reference_urls = [
            # "ftp://ftp.ncbi.nlm.nih.gov//genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz", 
            # "ftp://ftp.ncbi.nlm.nih.gov//genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz"
            "https://raw.githubusercontent.com/genome-sci/python_bioinfo_2018/master/1-1/reference/s288c.fna", 
            "https://raw.githubusercontent.com/genome-sci/python_bioinfo_2018/master/1-1/reference/s288c.gff"
        ]
        filenames = ["s288c.fna", "s288c.gff"]
        for reference_url, filename in zip(reference_urls, filenames):
            bs = os.path.basename(reference_url)
            save_path = os.path.join(self.reference_dir, bs)

            if (not os.path.exists(save_path)) or force_download:
                print(f"Downloading {save_path}")
                urllib.request.urlretrieve(reference_url, save_path)

                ## For ncbi data
                # target_file = os.path.join(self.reference_dir, filename)
                # with gzip.open(save_path, mode="rb") as gzip_file:
                #     with open(target_file, mode="wb") as decompressed_file:
                #         shutil.copyfileobj(gzip_file, decompressed_file)


    def dump(self, sra_file):
        cmd = [
            self.fastq_dump, "--split-files", 
            "--outdir", self.fastq_dir, sra_file
        ]
        subprocess.run(cmd)

        bs = os.path.splitext(os.path.basename(sra_file))[0]
        fastq1 = os.path.join(self.fastq_dir, f"{bs}_1.fastq")
        cmd = ["gzip", fastq1]
        subprocess.run(cmd)

        fastq2 = os.path.join(self.fastq_dir, f"{bs}_2.fastq")
        cmd = ["gzip", fastq2]
        subprocess.run(cmd)


    def first_fastqc(self):
        for i in range(453566, 453572):
            fastq1 = os.path.join(self.fastq_dir, f"SRR{i}_1.fastq.gz")
            cmd = [
                self.fastqc, "--nogroup", "-o", self.fastqc_dir_1, fastq1
            ]
            subprocess.run(cmd)

            fastq2 = os.path.join(self.fastq_dir, f"SRR{i}_2.fastq.gz")
            cmd = [
                self.fastqc, "--nogroup", "-o", self.fastqc_dir_1, fastq2
            ]
            subprocess.run(cmd)


    def trimming(self):
        trim_log_dir = os.path.join(self.trimming_dir, "log")
        os.makedirs(trim_log_dir, exist_ok=True)
        for i in range(453566, 453572):
            fastq1 = os.path.join(self.fastq_dir, f"SRR{i}_1.fastq.gz")
            fastq2 = os.path.join(self.fastq_dir, f"SRR{i}_2.fastq.gz")
            paired1 = os.path.join(self.trimming_dir, f"paired_SRR{i}_1.trim.fastq.gz")
            unpaired1 = os.path.join(self.trimming_dir, f"unpaired_SRR{i}_1.trim.fastq.gz")
            paired2 = os.path.join(self.trimming_dir, f"paired_SRR{i}_2.trim.fastq.gz")
            unpaired2 = os.path.join(self.trimming_dir, f"unpaired_SRR{i}_2.trim.fastq.gz")
            truseq3_pe_2 = os.path.join(self.trim_adapt_dir, "TruSeq3-PE-2.fa")
            log_file = os.path.join(trim_log_dir, f"log_SRR{i}.txt")
            
            cmd = [
                "java", "-jar", "-Xmx512m", self.trimmomatic, "PE", 
                "-threads", self.threads, "-phred33", 
                "-trimlog", log_file, 
                fastq1, fastq2, paired1, unpaired1, 
                paired2, unpaired2, 
                f"ILLUMINACLIP:{truseq3_pe_2}:2:30:10", 
                "LEADING:20", "TRAILING:20", "SLIDINGWINDOW:4:15", 
                "MINLEN:36"
            ]
            subprocess.run(cmd)


    def second_fastqc(self):
        for i in range(453566, 453572):
            fastq1 = os.path.join(self.trimming_dir, f"paired_SRR{i}_1.trim.fastq.gz")
            cmd = [
                self.fastqc, "--nogroup", "-o", self.fastqc_dir_2, fastq1
            ]
            subprocess.run(cmd)

            fastq2 = os.path.join(self.trimming_dir, f"paired_SRR{i}_2.trim.fastq.gz")
            cmd = [
                self.fastqc, "--nogroup", "-o", self.fastqc_dir_2, fastq2
            ]
            subprocess.run(cmd)


    def mapping(self):
        output_file = os.path.join(self.hisat_dir, os.path.basename(self.reference_fna))
        cmd = [
            self.hisat2_build, self.reference_fna, output_file
        ]
        subprocess.run(cmd)

        for i in range(453566, 453572):
            fastq1 = os.path.join(self.trimming_dir, 
                                  f"paired_SRR{i}_1.trim.fastq.gz")
            fastq2 = os.path.join(self.trimming_dir, 
                                  f"paired_SRR{i}_2.trim.fastq.gz")
            samfile = os.path.join(self.hisat_dir, f"SRR{i}.sam")
            bamfile = os.path.join(self.hisat_dir, f"SRR{i}.bam")
            sorted_bamfile = os.path.join(self.hisat_dir, f"SRR{i}.sorted.bam")
            cmd = [
                self.hisat2, "-p", self.threads, "-x", output_file, 
                "-1", fastq1, "-2", fastq2, "-S", samfile
            ]
            subprocess.run(cmd)

            cmd = [
                self.samtools, "view", "-@", self.threads, 
                "-b", samfile, ">", bamfile
            ]
            cmd = " ".join(cmd)
            subprocess.run(cmd, shell=True)

            cmd = [
                self.samtools, "sort", "-@", self.threads, 
                bamfile, ">", sorted_bamfile
            ]
            cmd = " ".join(cmd)
            subprocess.run(cmd, shell=True)


    def read_count(self):
        bs = os.path.splitext(os.path.basename(self.reference_gff))[0]
        reference_gene_id = os.path.join(
            self.reference_dir, f"{bs}_e.gff"
        )

        # add gene id
        records = list(GFF.parse(open(self.reference_gff)))
        gene_cnt = 0
        for r in records:
            for f in r.features:
                if f.type == "gene" or f.type == "pseudogene":
                    gene_cnt += 1
                    gene_id = "gene_" + str(gene_cnt).zfill(4)
                    f.qualifiers["gene_id"] = [gene_id]
                    for sf in f.sub_features:
                        sf.qualifiers["gene_id"] = [gene_id]
                        for ssf in sf.sub_features:
                            ssf.qualifiers["gene_id"] = [gene_id]

        with open(reference_gene_id, "w") as fh:
            GFF.write(records, fh)

        # Feature counts
        output_file = os.path.join(os.getcwd(), self.featurecounts_dir, "counts.txt")
        bam_files = [os.path.join(self.hisat_dir, f"SRR{i}.sorted.bam") for i in range(453566, 453572)]
        cmd = [
            self.featureCounts, "-p", "-T", self.threads, 
            "-t", "exon", "-g", "gene_id", 
            "-a", reference_gene_id, 
            "-o", output_file
        ] + bam_files
        subprocess.run(cmd)


    def setup(self, force_download=False):
        print("start setup software")
        self.setup_software(force_download)
        print("finish setup software")
        
        print()
        
        print("start download data")
        self.download_data(force_download)
        print("finish download")


    def run(self):
        print("start fastqc before trimming")
        self.first_fastqc()
        print("finish fastqc")

        print()

        print("start trimming")
        self.trimming()
        print("finish trimming")

        print()

        print("start fastqc after trimming")
        self.second_fastqc()
        print("finish fastqc")

        print()

        print("start mapping")
        self.mapping()
        print("finish mapping")

        print("start count")
        self.read_count()
        print("finish count")

        print("finish all process")


def main():
    fastq_dir = "sra"
    reference_dir = "reference"
    fastqc_dir_1 = "fastqc_before"
    fastqc_dir_2 = "fastqc_after"
    trimming_dir = "trimmomatic"
    hisat_dir = "hisat"
    counts_dir = "featurecount"
    tool_dir = "tools"
    fastq_dump = os.path.join(tool_dir, "sratoolkit.3.0.1-ubuntu64/bin/fastq-dump")
    trimmomatic = os.path.join(tool_dir, "Trimmomatic-0.39/trimmomatic-0.39.jar")
    fastqc = os.path.join(tool_dir, "FastQC/fastqc")
    hisat2_build = os.path.join(tool_dir, "tools/hisat2/hisat2-build")
    featureCounts = os.path.join(tool_dir, "subread-2.0.2-Linux-x86_64/bin/featureCounts")
    threads = 4

    sc = Scrna_seq(fastq_dir, reference_dir, 
                   fastqc_dir_1, fastqc_dir_2, 
                   trimming_dir, hisat_dir, 
                   counts_dir, tool_dir, 
                   threads=threads)
    
    mode = sys.argv[1]
    if mode == "setup":
        sc.setup()
    elif mode == "run":
        sc.run()
    else:
        print(f"Unknown argument : {mode}")
        print("Argument must be 'setup' or 'run'")
    


if __name__ == "__main__":
    main()
