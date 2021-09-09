# orcanizator
generate ORCA input files from SMILES

**Global Variables**
---------------
- **ORCA_EXECUTABLE**

---

## <kbd>function</kbd> `orca_generator`

```python
orca_generator(smiles: str, nb: int) → None
```

Generate ORCA input file from SMILES. ORCA arguments are hardcoded into the function. 



**Arguments:**
  smiles {str} -- SMILES structure of the molecule.  
  nb {int} -- number. This will be the name of the input file. 


---

## <kbd>function</kbd> `orca_reader`

```python
orca_reader(smiles: str, nb: int) → dict
```

Read following properties of ORCA output file: 
1. HOMO energy 
2. LUMO energy 
3. Dipole moment 
4. Isotropic polarizability 
5. Total thermal energy, U 
6. Total enthalpy, H 
7. Entropy term, S 
8. Gibbs free energy, G 
9. Rotational constant, A 
10. Rotational constant, B 
11. Rotational constant, C 



**Arguments:**
  smiles {str} -- SMILES structure of the molecule.  
  nb {int} -- number. This is the name of the output file to be read. 



**Returns:**
  {dict} -- a dictionary of read properties. 


---

## <kbd>function</kbd> `orcanize`

```python
orcanize(smiles: str, nb: int) → dict
```

This function does the following steps: 
1. Generate ORCA input file based on molecule SMILES by calling orca_generator function. 
2. Call ORCA job with generated input file. 
3. Read calculated properties from generated output file and return them as dictionary. This is done by calling orca_reader function. 

SMILES of successfully processed molecules are saved into done.log file. The done.log also contains the number which is input and output files name. If there was an error in calculations, then the molecule's SMILES is saved into undone.log file. 



**Arguments:**
  smiles {str} -- [SMILES molecule structure]  nb {int} -- [Number. This is ORCA input and output file name] 



**Returns:**
  {dict} -- [Dictionary of properties] 


---

## <kbd>function</kbd> `orcanize_many`

```python
orcanize_many(molecules: List[str], save: bool = True) → DataFrame
```

Run orcanize function with a collection of SMILES. 



**Arguments:**
  molecules {Union[str]} -- collection of SMILES strings.  
  save {bool} -- save the pandas dataframe to csv(True, default) or  not(False); by convention it is saved everytime ORCA  successfully finishes calculation 



**Returns:**
  {pd.DataFrame} -- pandas dataframe with read properties 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._ ... and then modified
