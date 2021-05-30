namespace MAP4

module MHFP =

    open NumpyDotNet
    open System
    open GraphMolWrap
    open MAP4.Utils

    let molToRingAtomLists (mol: ROMol) =
        use mutable rings = new Int_Vect_Vect()
        RDKFuncs.symmetrizeSSSR (mol, rings) |> ignore
        let list_of_rings = ivvToListList rings
        list_of_rings

    let isBondInMol (mol: ROMol) pair =
        let bond =
            mol.getBondBetweenAtoms (fst pair, snd pair)

        if bond <> null then
            Some(int (bond.getIdx ()))
        else
            None

    let bondsInRing (mol: ROMol) ring =
        let pairs =
            seq {
                for i in ring do
                    for j in ring do
                        if i <> j then uint32 i, uint32 j
            }

        pairs
        |> Seq.map (isBondInMol mol)
        |> Seq.choose id
        |> Set.ofSeq

    let shingleFromPath (mol: ROMol) bondIds =
        let subMol =
            RDKFuncs.pathToSubmol (mol, (listToVec bondIds))

        RDKFuncs.MolToSmiles(subMol)

    let shingleFromRing (mol: ROMol) ring =
        ring |> (bondsInRing mol) |> (shingleFromPath mol)

    let allAtomSmiles (mol: ROMol) =
        atomsSeq mol
        |> Seq.map RDKFuncs.GetAtomSmiles
        |> List.ofSeq

    let smilesFromEnvironment (mol: ROMol) (index: int) (environment: Int_Vect) =
        use mutable amap = new Int_Int_Map()

        try
            let submol =
                RDKFuncs.pathToSubmol (mol, environment, false, amap)

            if amap.Keys.Contains(index) then
                let smiles =
                    RDKFuncs.MolToSmiles(submol, true, false, amap.[index], true, false, false)

                if smiles.Length > 0 then
                    Some smiles
                else
                    None
            else
                None
        with :? ApplicationException -> None

    let singleAtomIdxEnvironmentShingle (mol: ROMol) (maxRadius: int) (atomIdx: int) =
        seq { 1 .. maxRadius }
        |> Seq.map (fun r -> RDKFuncs.findAtomEnvironmentOfRadiusN (mol, uint32 r, uint32 atomIdx))
        |> Seq.map (smilesFromEnvironment mol atomIdx)
        |> Seq.choose id
        |> Set.ofSeq

    let environmentShingles (mol: RWMol) (maxRadius: int) =
        seq { 0 .. int (mol.getNumAtoms ()) - 1 }
        |> Seq.map (singleAtomIdxEnvironmentShingle mol maxRadius)
        |> Seq.reduce (+)
        |> Set.toList

    let shinglingFromMol (mol: RWMol) (radius: int) (rings: bool) (minRadius: int) =

        let ringShingles =
            if rings then
                (molToRingAtomLists mol)
                |> List.map (shingleFromRing mol)
            else
                List.empty

        let singleAtomShingles =
            if minRadius = 0 then
                allAtomSmiles mol
            else
                List.empty

        let envShingles = environmentShingles mol radius
        ringShingles @ singleAtomShingles @ envShingles

    let shinglingFromMolDefaults mol = shinglingFromMol mol 3 true 1


    let shinglingFromSmiles (smiles: Smiles) (radius: int) (rings: bool) (minRadius: int) (sanitize: bool) =
        maybe {
            let! mol = smilesToMol smiles sanitize
            return shinglingFromMol mol radius rings minRadius
        }

    let shinglingFromSmilesDefaults smiles =
        shinglingFromSmiles smiles 3 true 1 true

    let hash (shingling: string list) =
        let hashValues = shingling |> List.map (stringToHash 4)
        np.array (List.toArray hashValues, dtype = np.UInt32)

    let fold (length: int) (hashValues: ndarray) =
        let mutable folded =
            np.zeros (shape (length), dtype = np.UInt8)

        let bits = np.``mod`` (hashValues, length)
        np.put (folded, bits, 1) |> ignore
        folded

    let secfpFromMol (mol: RWMol, length: int, radius: int, rings: bool, minRadius: int) =
        shinglingFromMol mol radius rings minRadius
        |> hash
        |> fold length

    let secfpFromSmiles (smiles: Smiles) (length: int) (radius: int) (rings: bool) (minRadius: int) (sanitize: bool) =
        maybe {
            let! mol = smilesToMol smiles sanitize
            return secfpFromMol (mol, length, radius, rings, minRadius)
        }
        
    let secfpFromSmilesDefault smiles =
        secfpFromSmiles smiles 2048 3 true 1 false

    let distance (a: ndarray) (b: ndarray) =
        let length = int a.shape.iDims.[0]
        let intersect = np.sum(np.equal (a, b)).item () :?> int
        1.0 - float (intersect) / float (length)

    type MHFP(?permutations, ?seed) =
        let permutations0 = defaultArg permutations 2048
        let seed0 = defaultArg seed 42
        let rand0 = np.random ()

        let mutable pa =
            np.zeros (shape (permutations0), np.UInt32)

        let mutable pb =
            np.zeros (shape (permutations0), np.UInt32)

        let maxHash0 = (1UL <<< 32) - 1UL

        do
            rand0.seed (seed0) |> ignore

            for i in [ 0 .. permutations0 - 1 ] do
                let mutable a =
                    rand0.randint (1L, maxHash0, shape (1), dtype = np.UInt32)

                let mutable b =
                    rand0.randint (0L, maxHash0, shape (1), dtype = np.UInt32)

                while unbox<bool> (np.isin(a, pa).item ()) do
                    a <- rand0.randint (1L, maxHash0, shape (1), dtype = np.UInt32)

                while unbox<bool> (np.isin(b, pb).item ()) do
                    b <- rand0.randint (0L, maxHash0, shape (1), dtype = np.UInt32)

                pa.[i] <- a
                pb.[i] <- b


        member this.prime = (1L <<< 61) - 1L
        member this.maxHash = maxHash0
        member this.permutations = permutations0
        member val seed = defaultArg seed 42
        member this.rand = rand0

        member this.permutations_a =
            pa.reshape (shape (this.permutations, 1))

        member this.permutations_b =
            pb.reshape (shape (this.permutations, 1))

        member this.permutations_a_64 =
            np.zeros (shape this.permutations, dtype = np.UInt64)

        member this.permutations_b_64 =
            np.zeros (shape this.permutations, dtype = np.UInt64)


        member this.fromMolecularShingling(tokens: Smiles list) =
            let mutable hashValues =
                np.full (shape (this.permutations, 1), this.maxHash, np.UInt32)

            for t in tokens do
                let th = stringToHash 4 t

                let hashes =
                    np.remainder (
                        np.remainder (this.permutations_a * th + this.permutations_b, this.prime),
                        this.maxHash
                    )

                hashValues <- np.minimum (hashValues, hashes)

            hashValues.reshape(shape (1, this.permutations)).[0] :?> ndarray

        member this.encode(smiles: Smiles) =
            maybe {
                let! tokens = shinglingFromSmiles smiles 3 true 1 true
                let hash = this.fromMolecularShingling (tokens)
                return hash
            }
