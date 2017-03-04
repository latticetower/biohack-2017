import pickle 

with open('set_pdb.pickle', 'rb') as f:
	set_pdb = pickle.load(f)
	with open("list_pdb.txt","w") as f2:
		for word in set_pdb:
			f2.write(word + ' ')