namespace MAP4

open System.Collections.Generic
open System
open System.Security.Cryptography
open GraphMolWrap
open NumpyDotNet
open FSharp.Collections
module Utils =
    type Smiles = string
    type FP =
    | Folded of ndarray
    | Unfolded of ndarray
        member this.getValue() = match this with
                                       | Folded a -> a
                                       | Unfolded a -> a

    let sha1 = HashAlgorithm.Create("SHA1")

    let encoding = System.Text.UTF8Encoding()

    let listToVec (l:int Set) =
        let result = new Int_Vect()
        l |> Set.iter result.Add
        result
        
    let intToVec(i: int) =
        let result = new Int_Vect()
        result.Add(i)
        result

    /// Create a seq of all atoms from a molecule
    let atomsSeq (mol:ROMol) =
        seq{for atom in mol.getAtoms() do yield atom}

    /// convert a string into a byte array using its utf-8 representation
    let stringToBytes (s: string) =
        encoding.GetBytes s
        
    /// Convert a byte array into a uint32    
    let bytesToInt (b:byte[]) endIndex =
        BitConverter.ToUInt32(b.[0..endIndex], 0)
    /// Use the first N bytes of the sha1 hash of a string and interpret them as uint32    
    let stringToHash (endAtByte:int) (s: string) =
        let bytes = stringToBytes s
        let hash = sha1.ComputeHash(bytes)
        bytesToInt hash endAtByte
          
    let ivvToListList (ivv: Int_Vect_Vect) =
        List.ofSeq(seq{ for i in ivv do List.ofSeq(seq{for j in i do j })})

    let randomArray(size:int, rand: np.random, minVal:int64, maxVal: uint64) =
        let h = HashSet()
        while h.Count < size do
            let v = rand.randint(minVal, Nullable(maxVal), newshape=shape(1), dtype=np.UInt32)
            h.Add(uint32 v) |> ignore
        let rand_list = Seq.fold (fun l se -> uint32(se)::l) [] h
        let result = np.array(Array.ofList(rand_list), dtype=np.UInt32)
        result.reshape([|size|])
        
    let smilesToMol(smiles: Smiles) (sanitize: bool) =
        try
            Some (RWMol.MolFromSmiles(smiles, 0, sanitize))
        with
        | e -> None
        
    /// Calculate the distance matrix for a molecule using the Floyd Warshall algorithm
    /// This is a simpler version compared to the RDKit implementation always using the edge weight 1.0 
    let distanceMatrix (mol:RWMol) =  // Floyd Warshall
        let numAtoms = int(mol.getNumAtoms())
        let dist = Array2D.create numAtoms numAtoms 10000000
        let numBonds = int(mol.getNumBonds())

        for i in [0..numBonds-1] do
            let bond = mol.getBondWithIdx(uint32(i))
            let fst = int(bond.getBeginAtomIdx())
            let snd = int(bond.getEndAtomIdx())
            dist.[fst,snd] <- 1
            dist.[snd,fst] <- 1
            
        let next = Array2D.create numAtoms numAtoms 0
        
        for k in [0..numBonds-1] do
            for i in [0..numBonds-1] do
                for j in [0..numBonds-1] do
                    if dist.[i,k] + dist.[k,j] < dist.[i,j] then
                        dist.[i,j] <- dist.[i,k] + dist.[k,j]
                        next.[i,j] <- next.[i,k]
        for i in [0..numAtoms-1] do
            dist.[i,i] <- 0
        dist       
        
    type MaybeBuilder() =
        member _.Bind(m, f) =
            Option.bind f m
            
        member _.Return(x) =
            Some x
            
        member _.ReturnFrom(x) =
            x
            
        member _.Zero() =
            None
            
    let rec comb n l = 
        match n, l with
        | 0, _ -> [[]]
        | _, [] -> []
        | k, (x::xs) -> List.map ((@) [x]) (comb (k-1) xs) @ comb k xs

    let maybe = MaybeBuilder()
        
            
        