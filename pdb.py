import pickle 
import requests



r = requests.get("http://www.genome.jp/dbget-bin/www_bget?hsa:10395")
page = r.text

count = 0
word = '1111' 
ok = 0
with open('set_pdb.pickle', 'rb') as f:
	set_pdb = pickle.load(f)
	for s in page:
		
		if (ok != 0):
			if (ok < 5):
				ok += 1
			else:
				set_pdb.add(word.upper())
				ok = 0
		word = word[1:]
		word += s
		l = len(word)
		if (word == 'pdb:'):
			ok = 1

	with open('set_pdb.pickle', 'wb') as f:
		pickle.dump(set_pdb, f)

