module MAP4.LSHForest

open System.Collections.Generic
open MAP4.Utils
open NumpyDotNet

type LSHForest(?dims: int, ?nPrefixTrees) =
    let mutable count = 0
    member val dims = defaultArg dims 2048
    member val nPrefixTrees = defaultArg nPrefixTrees 64
    member this.maxDepth = int(this.dims / this.nPrefixTrees)
    member this.hashTables = List.ofSeq (seq {for _ in [0..this.nPrefixTrees-1] do yield Dictionary<string, int>()})
    member this.hashRanges = np.array(Array.ofList(List.concat(seq{for i in [0..this.nPrefixTrees-1] do yield [i * this.maxDepth; (i+1) * this.maxDepth]}))).reshape(-1,2)
    member this.keys = Dictionary<int,int>()
    member this.clean = false
    
    
    member private this.swap(hashes:ndarray) =
        hashes.byteswap().rawdata
    
    member this.add(key: int , mhfp: ndarray) =
        this.keys.Add(0,0)
        


    
    
    
