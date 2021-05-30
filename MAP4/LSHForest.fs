module MAP4.LSHForest

open System
open System.Collections.Generic
open MAP4.Utils
open NumpyDotNet

type LSHForest(?dims: int, ?nPrefixTrees) =
    let mutable _count = 0
    let mutable _clean = false
    let dims = defaultArg dims 2048
    let nPrefixTrees = defaultArg nPrefixTrees 64
    let maxDepth = int(dims / nPrefixTrees)
    let hashTables = Array.ofSeq(seq {for _ in [0..nPrefixTrees-1] do yield Dictionary<ndarray, List<int>>()})
    let hashRanges = Array.ofSeq(seq{for i in [0..nPrefixTrees-1] do yield (i * maxDepth, (i+1) * maxDepth)})
    let keys = Dictionary<int,int>()
    member this.clean
        with get() = _clean
        and set(v) = _clean <- v
        
    member this.count
        with get() = _count
        and set(v) = _count <- v

    
    
    member private this.swap(hashes:ndarray) =
        hashes.byteswap()
    
    member this.add(key: int , mhfp: ndarray) =
        keys.Add(key, 0)
        let k = seq{for start, last in hashRanges do yield this.swap(mhfp.[Slice(start,last)] :?> ndarray) }
        for h, hashtable in Seq.zip k hashTables  do
            if not (hashtable.ContainsKey(h)) then
                hashtable.Add(h, List<int>())
            hashtable.[h].Add(key)
        this.clean <- false
            
        
        


    
    
    
