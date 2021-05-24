namespace MAP4.MHFP

open NumpyDotNet
open System
open GraphMolWrap
open MAP4.Utils

type MHFP(?permutations, ?seed) =
    let permutations0 = defaultArg permutations 2048
    let seed0 = defaultArg seed 42
    let rand0 = np.random()
    
    let mutable pa = np.zeros(shape(permutations0), np.UInt32)
    let mutable pb = np.zeros(shape(permutations0), np.UInt32)
    let maxHash0 = (1UL <<< 32) - 1UL
    
    do
        rand0.seed(seed0) |> ignore
        for i in ([0..permutations0-1]) do
            let mutable a = rand0.randint(1L, maxHash0, shape(1), dtype=np.UInt32)
            let mutable b = rand0.randint(0L, maxHash0, shape(1), dtype=np.UInt32)
            while unbox<bool>(np.isin(a, pa).item()) do
                a <- rand0.randint(1L, maxHash0, shape(1), dtype=np.UInt32)
            while unbox<bool>(np.isin(b, pb).item()) do
                b <- rand0.randint(0L, maxHash0, shape(1), dtype=np.UInt32)
            pa.[i] <- a
            pb.[i] <- b
               
    
    member this.prime = (1L <<< 61) - 1L
    member this.maxHash = maxHash0
    member this.permutations = permutations0
    member val seed = defaultArg seed 42
    member this.rand = rand0

    member this.permutations_a = pa.reshape(shape(this.permutations, 1))
    member this.permutations_b = pb.reshape(shape(this.permutations, 1))         
    member this.permutations_a_64 = np.zeros (shape this.permutations, dtype = np.UInt64)
    member this.permutations_b_64 = np.zeros (shape this.permutations, dtype = np.UInt64)

    

    static member molToRingAtomLists(mol: ROMol) =
        use mutable rings = new Int_Vect_Vect()
        RDKFuncs.symmetrizeSSSR (mol, rings) |> ignore
        let list_of_rings = ivvToListList rings
        list_of_rings

    static member isBondInMol (mol: ROMol) pair =
        let bond = mol.getBondBetweenAtoms (fst pair, snd pair)
        if bond <> null then
            Some(int (bond.getIdx ()))
        else
            None

    static member bondsInRing (mol: ROMol) ring =
        let pairs =
            seq {for i in ring do for j in ring do if i <> j then uint32 i, uint32 j}
        pairs
        |> Seq.map (MHFP.isBondInMol mol)
        |> Seq.choose id
        |> Set.ofSeq

    static member shingleFromPath (mol: ROMol) bondIds =
        let subMol =
            RDKFuncs.pathToSubmol (mol, (listToVec bondIds))

        RDKFuncs.MolToSmiles(subMol)

    static member shingleFromRing (mol: ROMol) ring =
        ring
        |> (MHFP.bondsInRing mol)
        |> (MHFP.shingleFromPath mol)

    static member allAtomSmiles(mol: ROMol) =
        atomsSeq mol
        |> Seq.map RDKFuncs.GetAtomSmiles
        |> List.ofSeq

    static member smilesFromEnvironment (mol: ROMol) (index: int) (environment: Int_Vect) =
        use mutable amap = new Int_Int_Map()
        try
            let submol = RDKFuncs.pathToSubmol (mol, environment, false, amap)
            if amap.Keys.Contains(index) then
                let smiles = RDKFuncs.MolToSmiles(submol, true, false, amap.[index], true, false, false)
                if smiles.Length > 0 then Some smiles else None
            else
                None
        with :? ApplicationException -> None

    static member singleAtomIdxEnvironmentShingle (mol: ROMol) (maxRadius: int) (atomIdx: int) =
        seq { 1 .. maxRadius }
        |> Seq.map (fun r -> RDKFuncs.findAtomEnvironmentOfRadiusN (mol, uint32 r, uint32 atomIdx))
        |> Seq.map (MHFP.smilesFromEnvironment mol atomIdx)
        |> Seq.choose id
        |> Set.ofSeq

    static member environmentShingles (mol: RWMol) (maxRadius: int) =
        seq { 0 .. int (mol.getNumAtoms()) - 1 }
        |> Seq.map (MHFP.singleAtomIdxEnvironmentShingle mol maxRadius)
        |> Seq.reduce (+)
        |> Set.toList

    static member shinglingFromMol(mol: RWMol, ?radius: int, ?rings: bool, ?minRadius: int) =
        let radius = defaultArg radius 3
        let rings = defaultArg rings true
        let minRadius = defaultArg minRadius 1

        let ringShingles =
            if rings then
                (MHFP.molToRingAtomLists mol)
                |> List.map (MHFP.shingleFromRing mol)
            else
                List.empty

        let singleAtomShingles =
            if minRadius = 0 then
                MHFP.allAtomSmiles mol
            else
                List.empty

        let envShingles = MHFP.environmentShingles mol radius
        ringShingles @ singleAtomShingles @ envShingles
        
        
    static member shinglingFromSmiles(smiles: Smiles, ?radius: int, ?rings: bool, ?minRadius: int, ?sanitize: bool) =
        let sanitize = defaultArg sanitize true
        let minRadius = defaultArg minRadius 1
        let rings = defaultArg rings true
        let radius = defaultArg radius 3

        maybe {
            let! mol = smilesToMol smiles sanitize
            return MHFP.shinglingFromMol (mol, radius, rings, minRadius)
        }

    member this.fromMolecularShingling (tokens: Smiles list) =
        let mutable hashValues = np.full(shape(this.permutations, 1), this.maxHash, np.UInt32)
       
        for t in tokens do
            let th = stringToHash 4 t
            let hashes = np.remainder (
                              np.remainder (
                                this.permutations_a * th + this.permutations_b, this.prime),
                              this.maxHash)
            hashValues <- np.minimum(hashValues, hashes)

        hashValues.reshape(shape (1, this.permutations)).[0] :?> ndarray

    member this.encode(smiles: Smiles) =
        maybe {
            let! tokens = MHFP.shinglingFromSmiles smiles
            let hash = this.fromMolecularShingling tokens
            return hash
        }
        
    static member hash(shingling: string list) =
        let hashValues = shingling |> List.map (stringToHash 4)
        np.array(List.toArray hashValues, dtype=np.UInt32)
        
    static member fold (length: int) (hashValues: ndarray) =
        let mutable folded = np.zeros(shape(length), dtype=np.UInt8)
        let bits = np.``mod``(hashValues, length)
        np.put(folded, bits,1) |> ignore
        folded
        
    static member secfpFromMol(mol: RWMol, ?length : int, ?radius:int, ?rings: bool, ?minRadius: int) =
        let length = defaultArg length 2048
        let radius = defaultArg radius 3
        let rings = defaultArg rings true
        let minRadius = defaultArg minRadius 1
        MHFP.shinglingFromMol(mol, radius, rings, minRadius) |> MHFP.hash |> MHFP.fold length
        
    static member secfpFromSmiles(smiles: Smiles, ?length : int, ?radius:int, ?rings: bool, ?minRadius: int,
                                  ?sanitize: bool) =
        let length = defaultArg length 2048
        let radius = defaultArg radius 3
        let rings = defaultArg rings true
        let minRadius = defaultArg minRadius 1
        let sanitize = defaultArg sanitize false
        
        maybe {
            let! mol = smilesToMol smiles sanitize
            return MHFP.secfpFromMol(mol, length, radius, rings, minRadius)
        }
        
    static member distance (a : ndarray) (b: ndarray) =
        let length = int a.shape.iDims.[0]
        let intersect = np.sum(np.equal(a, b)).item() :?> int
        1.0 - float(intersect) / float(length)


        
