import urllib.request

url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.txt"
output_file = "20130606_sample_info.txt"

urllib.request.urlretrieve(url, output_file)
print("Download complete.")
