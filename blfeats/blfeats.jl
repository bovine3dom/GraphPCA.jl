module blFeats

    import LightGraphs
    import JLD2
    import FileIO
    import YAML
    import DataFrames
    const df = DataFrames
    const lg = LightGraphs
    include("../offline-laptop/mykavosh/kavosh.jl")

    macro p(s)
        quote
            $(esc(parse("fetch.([@spawn " * "$s"[2:end] * ")")))
        end
    end

    # Todo
    # - sensible names / categories
    # - sensible directory structure
    # - sensible way of loading / generating from source files
    # - convert graphholder etc. to new system
    # - datamat should accept a list of features
    # - should parallelise only if feature doesn't exist?
    # - make all keys case insensitive
    # - should probably use "get(dict,key,default)" instead of accessing dict directly

    const CACHE_DIR = "/home/olie/.dropbox-work/julia_graph_cache/"

    type Feature
        name::String
        padfunction::Function # one arg, graph und meta[]
        getfunction::Function # with one argument, graph und meta
    end

    FEATURE_DICT = Dict()

    type GraphAndMeta
        graph::LightGraphs.SimpleGraphs.AbstractSimpleGraph
        meta::Dict{String,Any}
        feature::Dict{String,Any}
    end

    function GraphAndMeta(graph; name=randstring(32))
        meta = Dict(
            "name" => name,
            "location" => CACHE_DIR,
            "category" => "Uncategorised",
        )
        return GraphAndMeta(graph,meta,Dict())
    end

    function konect_loader(directory)
        meta_file = split(readstring(`fish -c "ls $(directory)/**meta.*"`))[1]
        graph_file = split(readstring(`fish -c "ls $(directory)/**out.*"`))[1]
        meta = Dict()
        # Handy script for tidying up metadata:
        # sed 's/:/-/2' $(filename)
        # should probably do it in place here rather than messing with the source files
        
        try
            meta = YAML.load_file(meta_file)
            if haskey(meta,"Name") meta["name"] = meta["Name"] end
        catch(error)
            if isa(error,YAML.ParserError) || isa(error,YAML.ScannerError)
                meta["name"] = split(graph_file,"/")[end]
            else
                println(typeof(error))
                throw(error)
            end
        end
        meta["location"]=CACHE_DIR
        
        # Should change type of graph we make based on this
        meta["type"] = join(split(strip(readstring(`head -n1 $(graph_file)`)), " ")[2:end], " ")
        
        cachedfile = meta["location"]*meta["name"]*".jld2"
        try
            candidate = loadGraphAndMeta(cachedfile)
            merge!(candidate.meta,meta)
            return candidate
        catch(e)
            if isa(e,SystemError)
                print(e)
                graph = graph_from_konect(graph_file)
                return GraphAndMeta(graph,meta,Dict())
            else
                throw(e)
            end
        end
    end

    function graph_from_konect(filename)
        edgeframe= df.readtable(filename,separator=' ',allowcomments=true,commentmark='%')
        edgelist = Array(edgeframe[:,1:2])
        numnodes = maximum(Array(edgelist[:,1:2])[:])
        a = lg.Graph(numnodes)
        for i in 1:size(edgelist)[1]
            edge = tuple(edgelist[i,1:2]...)
            lg.add_edge!(a,edge)
        end
        return a
    end

    function thomps_loader(filename,directed=true)
        meta = Dict()
        meta["location"]=CACHE_DIR
        meta["name"] = split(split(filename,'/')[end],"txt.")[1]
        meta["name"] = directed ? meta["name"] : meta["name"]*"_undirected"
        meta["category"] = "Animal"
        meta["cite"] = "Thompson & Townsend"
        meta["type"] = "sym unweighted"
        cachedfile = meta["location"]*meta["name"]*".jld2"
        try
            candidate = loadGraphAndMeta(cachedfile)
            merge!(candidate.meta,meta)
            return candidate
        catch(e)
            if isa(e,SystemError)
                print(e)
                graph = graph_from_thomps(filename)
                graph = directed ? graph : lg.Graph(graph)
                return GraphAndMeta(graph,meta,Dict())
            else
                throw(e)
            end
        end
    end

    function graph_from_thomps(filename)
        raw = df.readtable(filename,separator='\t')
        adj = Matrix(raw[:,2:end])
        #adj = adj'*adj # not a real undirected graph, want binary or with transpose
        lg.DiGraph(adj)
    end

    #konect_graphs = []
    #konect_strings = split(readstring(`fish -c "ls /home/olie/Dropbox_work/big_data/data_oct16/konect/**/out.* "`))

    function generategraph(graphfunc,args)
        graph = GraphAndMeta(graphfunc(args...),name=string(graphfunc)*join(string.(args),"_")*"_"*randstring(32))
        graph.meta["function"] = graphfunc
        graph.meta["args"] = args
        graph.meta["category"] = "Generated"
        return graph
    end
         

    function getfeature!(gm,feature::String;)
        if haskey(gm.feature,feature) return gm.feature[feature] end
        gm.feature[feature] = FEATURE_DICT[feature].getfunction(gm)
        saveGraphAndMeta(gm)
        return gm.feature[feature]
    end

    function saveGraphAndMeta(gm;location="")
        location = location == "" ? gm.meta["location"]*gm.meta["name"]*".jld2" : location
        FileIO.save(location,gm.meta["name"],gm)
    end
        
    function loadGraphAndMeta(location)
        return [v for v in values(FileIO.load(location))][1]
    end

    function loadbyname(name)
        location = CACHE_DIR*name
        return loadGraphAndMeta(location)
    end

    type MotifFeature
        deets::Dict{Array{UInt64,1},Float64}
    end

    type AdjacencyFeature
        deets::SparseMatrixCSC{Int64,Int64}
    end

    function array_from_motif_dicts(motif_dicts)
        function mapping_from_motif_count_array(motifs)
            karr = []
            for ks in keys.(motifs)
                karr = vcat(karr, collect(ks))
            end
            #nauty.label_to_adj.(collect(Set(karr)),4)
            dmap = Dict()
            [dmap[j]=i for (i,j) in enumerate(collect(Set(karr)))]
            return dmap
        end

        dmap = mapping_from_motif_count_array(motif_dicts);
        result = Array{Float64,2}(length(dmap),0)
        for example in motif_dicts
            temparr = zeros(Float64,length(dmap),1)
            for (key, value) in example
                temparr[dmap[key]] = value
            end
            result = hcat(result,temparr)
        end
        rmap = Dict()
        for (k,v) in dmap
            rmap[v] = k
        end
        (result',dmap,rmap)
    end

    function getdatamat(graphs,feature::String)
        return FEATURE_DICT[feature].padfunction(graphs)
    end

    function loadall(filter="*")
        loadGraphAndMeta.(split(readstring(`fish -c "ls \"$(CACHE_DIR)\"$(filter)"`),"\n")[1:end-1])
    end

    function dictcat(dicts)
        ansdict = Dict()
        for dict in dicts
            for (k,v) in dict
                ansdict[k] = push!(get(ansdict,k,[]),v)
            end
        end
        ansdict
    end

    function dictstat(dicts)
        ansdict = Dict()
        catdict = dictcat(dicts)
        dictdot(catdict,v->(mean(v),std(Float64.(v))))
    end


    function zscore(x,μ,σ)
        (x-μ)/σ
    end

    function dictdot!(dict,func)
        for (k,v) in dict
            dict[k] = func(v)
        end
        dict
    end

    function dictdot(dict,func)
        ansdict = Dict()
        for (k,v) in dict
            ansdict[k] = func(v)
        end
        ansdict
    end

    function zmotif(motifdict,normdict)
        ansdict = Dict()
        for (k,v) in motifdict
            # rule of succession, sort of:
            # if motif does not show up at all, assume 1% prevelance
            # with std of 1pp.
            # really we need number of motifs to do this properly
            meanandstd = get(normdict,k,(0.01,0.01))
            ansdict[k] = zscore(v,meanandstd[1],meanandstd[2])
        end
        ansdict
    end

    function getclosedense(g,graphs)
        tol = 0.0001
        firsttry = filter(x-> isapprox(lg.density(x.graph),lg.density(g.graph),atol=tol), graphs)
        while length(firsttry) < 10
            tol *= 2
            firsttry = filter(x-> isapprox(lg.density(x.graph),lg.density(g.graph),atol=tol), graphs)
        end
        firsttry
    end

    function normmotifs(g,graphs)
        motifs = dictstat([blFeats.getfeature!(x,"motif_4") for x in getclosedense(g,graphs)])
        zmotif(blFeats.getfeature!(g,"motif_4"),motifs)
    end

    benchgraphs = loadall("*erdos*")

    function normmotifs(g)
        normmotifs(g, benchgraphs)
    end

    append_to_feature_dict!(feat::Feature) = FEATURE_DICT[feat.name] = feat

    append_to_feature_dict!(Feature(
            "motif_3",
            (graphs)->array_from_motif_dicts(
                [getfeature!(graph,"motif_3") for graph in graphs]),
            (g)->kavosh.getsubgraphs(g.graph,3))
    )

    append_to_feature_dict!(Feature(
            "motif_4",
            (graphs)->array_from_motif_dicts(
                [getfeature!(graph,"motif_4") for graph in graphs]),
            (g)->kavosh.getsubgraphs(g.graph,4))
    )
    
    append_to_feature_dict!(Feature(
            "motif_4_norm",
            (graphs)->array_from_motif_dicts(
                [getfeature!(graph,"motif_4_norm") for graph in graphs]),
            normmotifs)
    )

    append_to_feature_dict!(Feature(
            "motif_5",
            (graphs)->array_from_motif_dicts(
                [getfeature!(graph,"motif_5") for graph in graphs]),
            (g)->kavosh.getsubgraphs(g.graph,5))
    )

    append_to_feature_dict!(Feature(
            "motif_6",
            (graphs)->array_from_motif_dicts(
                [getfeature!(graph,"motif_6") for graph in graphs]),
            (g)->kavosh.getsubgraphs(g.graph,6))
    )
    # Todo: graph mixing matrices / sorted degree matrices
end
