def get_max_atg(orf_data, locus_model):

    max_orf_length = max(orf_data[1])
    max_orf_models = []
    max_orf_atg = []
    for i, orf_length in enumerate(orf_data[1]):
        if orf_length == max_orf_length:
            transcript_id = orf_data[0][i]
            max_orf_models.append(transcript_id)
            atg = locus_model.transcript_dict[transcript_id].get_atg()
            max_orf_atg.append(atg)
    assert len(set(max_orf_atg)) == 1
    return max_orf_atg[0], max_orf_models


def check_atg_universal(atg, locus_model):
    # Currently this method is disabled, by always returning True.
    transcript_dict = locus_model.transcript_dict
    for transcript_id in transcript_dict:
        atg_present = False
        transcript_model = transcript_dict[transcript_id]
        exon_list = transcript_model.exon_list
        for exon in exon_list:
            if atg in range(exon[0], exon[1]+1):
                atg_present = True
        if atg_present is True:
            return True

    return True


def characterise_max_orfs2(locus_dict):
    ambiguous_locus = []
    max_orf_transcripts = []
    max_orf_counts = []
    for locus_id in locus_dict:
        locus_model = locus_dict[locus_id]
        orf_data = locus_model.get_orf_lengths()

        # Some list contain None values, which is not valid for max() used ahead
        temp_lst = [n if n is not None else 0 for n in orf_data[1]]
        orf_data[1] = temp_lst
        
        atg, max_orf_models = get_max_atg(orf_data, locus_model)
        # Currently check_atg_universal is disabled, by always returning True.
        while check_atg_universal(atg, locus_model) is False:
            #remove all instances of current atg
            new_orf_data = [[], []]
            try:
                max_orf_length = max(orf_data[1])
            except: 
                break
            for i, orf_length in enumerate(orf_data[1]):
                if orf_length != max_orf_length:
                    new_orf_data[0].append(orf_data[0][i])
                    new_orf_data[1].append(orf_data[1][i])
            orf_data = new_orf_data
            try:
                atg, max_orf_models = get_max_atg(orf_data, locus_model)
            except:
                ambiguous_locus.append(locus_id)

        locus_model.rep_atg = atg
        locus_dict.update({locus_id: locus_model})

        max_orf_transcripts.extend(max_orf_models)
        n_max_models = len(max_orf_models)
        max_orf_counts.append(n_max_models)

    print("No universal ATG among annotated CDS models: ", len(ambiguous_locus))
        
    return max_orf_transcripts, max_orf_counts, locus_dict
