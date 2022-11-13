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
    let hashTables = Array.ofSeq(seq {for _ in [0..nPrefixTrees-1] do yield Dictionary<byte [], List<int>>()})
    let hashRanges = Array.ofSeq(seq{for i in [0..nPrefixTrees-1] do yield (i * maxDepth, (i+1) * maxDepth)})
    let hastablesSorted: ndarray array = Array.create nPrefixTrees (np.empty(1, np.Float32))
    let keys = Dictionary<int,int>()
    member this.clean
        with get() = _clean
        and set(v) = _clean <- v
        
    member this.count
        with get() = _count
        and set(v) = _count <- v
        
    member private this.swap(hashes:ndarray) =
        hashes.byteswap().tobytes()
    
    member this.add(key: int , mhfp: ndarray) =
        keys.Add(key, 0)
        let k = [|for start, last in hashRanges do yield this.swap(mhfp.[Slice(start,last)] :?> ndarray) |]
        for h, hashtable in Array.zip k hashTables  do
            if not (hashtable.ContainsKey(h)) then
                hashtable.Add(h, List<int>())
            hashtable.[h].Add(key)
        this.clean <- false
    
    member this.index() =
        for i,hashtable in Array.indexed hashTables do
            let a: obj[] = [|for k in hashtable.Keys do yield k|]
            let b = np.array(a)
            hastablesSorted.[i] <- b.Sort()
        this.clean <- true
        
    member this.binarySearch(n: int, hashtableSorted: ndarray, prefix: byte[], lenPrefix: int) =
        let mutable i = 0
        let mutable j = n
        while i < j do
            let h = int(i + (j - i) / 2)
            if ((hashtableSorted.[h] :?>ndarray).[Slice(lenPrefix)] :?> byte[]) >= prefix then
                i <- h + 1
            else
                j <- h
        
        
    member private this.internalQuery(mhfp: ndarray, r: int) =
        let prefixes = [ for start, _ in hashRanges do yield this.swap(mhfp.[Slice(start,start+r)] :?> ndarray) ]
        let lenPrefixes = prefixes.[0].Length
        1

        
        
        


    
    
    
