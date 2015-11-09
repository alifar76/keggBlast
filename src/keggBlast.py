import re
from urllib2 import urlopen
from ClientForm import ParseResponse
from bs4 import BeautifulSoup

from Bio import SeqIO
from datetime import datetime

def getProt(protein_seq):
	seq_data = open("sequence.gb", 'rU')
	records = SeqIO.parse(seq_data, "gb")
	for x in records:
		for y in x.features:
			try:
				protein_seq[y.qualifiers['protein_id'][0]] = str(y.qualifiers['translation'][0])
			except KeyError:
				pass
	return protein_seq

def keggBlast(proteins,filename):
	outfile = open(filename, "w")
	ko_prot = {}
	for protein in proteins.values():
		prot_id = [key for key, value in proteins.iteritems() if value == protein][0]
		kos = []
		response = urlopen('http://www.genome.jp/tools/blast/')
		forms = ParseResponse(response, backwards_compat=False)
		form = forms[0]
		form["sequence"] = protein
		filehand = urlopen(form.click()).read()
		soupy = BeautifulSoup(filehand)
		for link in soupy.find_all('a'):
			name_id = re.compile(r"http://www.genome.jp/dbget-bin/www_bget?").search(str(link))
			if name_id:
				url_name = str(link).split("<a href=")[1].split("</a>")[0].split('">')
				result_url = url_name[0].split('"')[1]
				filehandle = urlopen(result_url)
				soup = BeautifulSoup(filehandle)
				for link in soup.find_all('a'):
					ko_id = re.compile(r"ko\:(\K\d\d\d\d\d)").search(str(link))
					if ko_id:
						#print ko_id.group(1)+"\t"+url_name[1]+"\t"+prot_id
						#outfile.write(ko_id.group(1)+"\t"+url_name[1]+"\t"+prot_id+"\n")
						kos.append(ko_id.group(1))
		unique = list(set(kos))
		ko_prot[prot_id] = tuple(unique)
	for x in ko_prot.keys():
		for y in ko_prot[x]:
			outfile.write(x+"\t"+y+"\n")
	outfile.close()
	return outfile
					
if __name__ == '__main__':
	startTime = datetime.now()		
	proteins = getProt({})
	final_list = keggBlast(proteins,"koList.txt")
	print "Completion Time: "+ str(datetime.now()-startTime)	
