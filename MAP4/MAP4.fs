namespace MAP4.MAP4
open System
open System.Collections.Generic
open GraphMolWrap
open MAP4
open MAP4.MHFP
open Utils

type MAP4Calculator(?dimensions:int, ?radius:int, ?isCounted:bool, ?returnString:bool) =
    member _.dimensions = defaultArg dimensions 1024
    member _.radius = defaultArg radius 2
    member _.isCounted = defaultArg isCounted false
    member _.returnString = defaultArg returnString false
    member this.encoder = MHFP(this.dimensions)
    
    member private this.getAtomEnvs(mol:RWMol) =
        let ids = atomsSeq mol |> Seq.map (fun at -> int(at.getIdx()))
        let envs = ids |> Seq.map (singleAtomIdxEnvironmentShingle mol this.radius)
        Seq.zip ids envs |> Map.ofSeq
    
    member private this.allPairs(mol:RWMol, atomEnvs:Map<int,Set<string>>) =
        let atomPairs = List<string>()
        let distanceMatrix = distanceMatrix mol
        let numAtoms = int(mol.getNumAtoms())
        let shingleDict=Dictionary<string,int>()
        
        for ids in (comb 2 [0..numAtoms-1]) do
            let idx1 = ids.[0]
            let idx2 = ids.[1]
            let dist = int(distanceMatrix.[idx1, idx2]).ToString()

            for i in [0..this.radius-1] do
                let envA = (Set.toList atomEnvs.[idx1]).[i]
                let envB = (Set.toList atomEnvs.[idx2]).[i]
                let ordered = List.sort (envA :: envB :: [])
                let mutable shingle = $"{ordered.[0]}|{dist}|{ordered.[1]}"
                if this.isCounted then
                    if not(shingleDict.ContainsKey(shingle)) then
                        shingleDict.[shingle] <- 0
                    shingleDict.[shingle] <- shingleDict.[shingle] + 1 
                    shingle <- $"|{shingleDict.[shingle]}"
                atomPairs.Add shingle
        Set.toList(Set.ofArray(atomPairs.ToArray()))
        
    member private this.fold(pairs) =
        let fpHash = hash pairs
        Folded (fold this.dimensions fpHash)
        
    member this.calculateFolded(mol) =
        let atomEnvPairs = this.allPairs(mol, this.getAtomEnvs(mol))
        this.fold atomEnvPairs

    member this.calculate(mol) =
        let shingles = this.allPairs(mol, this.getAtomEnvs(mol))
        Unfolded (this.encoder.fromMolecularShingling shingles)
        
    member this.calculateMany(mols:RWMol list) =
        Array.Parallel.map this.calculate (List.toArray mols)
        
    member this.calculateManyFolded(mols: RWMol list) =
        Array.Parallel.map this.calculateFolded (List.toArray mols)
        
    member this.distance(fpA: FP, fpB: FP) =
        let dist =
            match fpA, fpB with
            | Folded a, Folded b -> MHFP.distance a b
            | Unfolded a, Unfolded b -> MHFP.distance a b
            | _, _ -> raise(ArgumentException("both arguments must be folded or unfolded fingerprints"))
        dist 

        

        
    

        
        
    

        
