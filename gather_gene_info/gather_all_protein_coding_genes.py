# this file will grab all the protein coding genes

import pandas as pd
import requests
from bs4 import BeautifulSoup
import random
import regex as re
from pickle import dump,load
import sys
import sqlite3
from collections import Counter
from os.path import exists
from multiprocessing import Pool

sys.setrecursionlimit(100000)

problem_file = "gathered_data/problem_genes.list"

def get_bp_and_cds(gene,debug=False):
  result = {}
  payload={}

  UAS = ("Mozilla/5.0 (Windows NT 6.1; WOW64; rv:40.0) Gecko/20100101 Firefox/40.1", 
          "Mozilla/5.0 (Windows NT 6.3; rv:36.0) Gecko/20100101 Firefox/36.0",
          "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10; rv:33.0) Gecko/20100101 Firefox/33.0",
          "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36",
          "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2227.1 Safari/537.36",
          "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2227.0 Safari/537.36",
          )

  ua = UAS[random.randrange(len(UAS))]
  headers = {
    'user-agent':ua
  }

  response1 = requests.request("GET", "https://www.ncbi.nlm.nih.gov/gene/?term={}".format(gene), headers=headers, data=payload)
  soup1 = BeautifulSoup(response1.text, 'html.parser')
  titles = soup1.findAll("a",id="feat_gene_title")
  result["gene"] = titles[0].contents[0].split()[0]
  taxons = soup1.findAll("p",class_="ncbi-doc-authors")
  result["organism"] = taxons[0].find("i").contents[0]
  result["gene_link"] = titles[0]["href"]
  
  response2 = requests.request("GET", titles[0]["href"], headers=headers, data=payload)
  soup2 = BeautifulSoup(response2.text, 'html.parser')
  result["gene_bank_url"] = "https://www.ncbi.nlm.nih.gov/"+soup2.findAll("a",title="Nucleotide GenBank report")[0]["href"]
  
  response3 = requests.request("GET", result["gene_bank_url"], headers=headers, data=payload)
  soup3 = BeautifulSoup(response3.text,"html.parser")
  soup3.findAll("meta",{"name":"ncbi_uidlist"})[0]["content"]
  result["ncbi_id"] = soup3.findAll("meta",{"name":"ncbi_uidlist"})[0]["content"]
  result["start"],result["stop"] = re.search(re.compile("from=(\d+)&to=(\d+)"),result["gene_bank_url"]).groups()
  result["strand"] = "on" if "strand=true" in result["gene_bank_url"] else "off"
  result["ncbi_phid"] = soup3.findAll("meta",{"name":"ncbi_phid"})[0]["content"].split()[0]
  
  target_url_template = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={}&db=nuccore&report=genbank&conwithfeat=on&basic_feat=on&hide-cdd=on&from={}&to={}&strand={}&retmode=html&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000"
  #target_url_template2 = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={}&db=nuccore&report=genbank&=&from={}&to={}&retmode=html&ncbi_phid={}&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000"
  result["genbank_jquery"] = target_url_template.format(result["ncbi_id"],result["start"],result["stop"],result["strand"])
  if debug:
    print(result["gene_bank_url"])
    print(result["genbank_jquery"])
  response4 = requests.request("GET", result["genbank_jquery"], headers=headers, data=payload)
  if debug:
    print(len(response4.text))
  m = re.search(re.compile(r"CDS\s+join\(((?:&|l|t|;|\.|,|\d|\s)+)\)"),response4.text)
  result["CDS"] = "".join(m.groups()[0].split())
  soup4 = BeautifulSoup(response4.text,"html.parser")
  result["seq"] = "".join(["".join(s.contents[0].split()) for s in soup4.findAll("span",class_="ff_line")])
  if debug:
    print("done!")
  return result

def get_all_human_protein_coding_genes():
  connection = sqlite3.connect("../data/EnsDb.Hsapiens.v79.sqlite")
  cursor = connection.cursor()
  rows = cursor.execute("SELECT gene_name,gene_biotype FROM gene").fetchall()
  protein_coding_genes = [g for g,b in rows if b == "protein_coding"]
  return protein_coding_genes

def get_all_gene_data(genes):
  for gene in genes:
    try:
      print(gene)
      result = get_bp_and_cds(gene)
      with open("gathered_data/raw_gene_results/{}.pickle","wb+".format(gene)) as f:
        dump(result,f)
    except Exception as e:
      print("failed")
      print(e)
      with open(problem_file,"a+") as f:
        f.write("{}\n".format(gene))

def get_gene_seq(gene):
    msg = ""
    if not exists("gathered_data/raw_gene_results/{}.pickle".format(gene)):
        msg += "Could not find data for {}\n".format(gene)
        try:
            result = get_bp_and_cds(gene)
            with open("gathered_data/raw_gene_results/{}.pickle".format(gene),"wb+") as f:
                dump(result,f)
            msg += "done!"
        except Exception as e:
            msg += "failed\n{}".format(e)
            with open(problem_file,"a+") as f:
                f.write("{}\n".format(gene))
    else:
        msg += "data for {} exists".format(gene)
    return msg

if __name__ == "__main__":
  connection = sqlite3.connect("../data/EnsDb.Hsapiens.v79.sqlite")
  cursor = connection.cursor()
  rows = cursor.execute("SELECT gene_name,gene_biotype FROM gene").fetchall()
  protein_coding_genes = [g for g,b in rows if b == "protein_coding"]
  with Pool(processes=32) as P:
    msgs = P.map(get_gene_seq, protein_coding_genes)
  with open("gathered_data/msgs.out","w+") as f:
    f.write("\n".join(msgs))
  obtained = sum("done" in m or "exists" in m for m in msgs)
  print("Obtained data on {}/{} ({:.2f}%) of protein-coding genes".format(obtained,len(protein_coding_genes),100*obtained/len(protein_coding_genes)))