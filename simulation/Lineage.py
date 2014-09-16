class NoAF_OutOfAfricaModel(MultiStageModel):
    '''A dempgrahic model for the CHB, CEU, and YRI populations, as defined in
    Gutenkunst 2009, Plos Genetics. The model is depicted in Figure 2, and the 
    default parameters are listed in Table 1 of this paper. The AF population is
    removed from the model after it splits into AF and B.'''
    def __init__(self, 
        T0,
        N_A=7300,
        N_AF=12300,
        N_B=2100,
        N_EU0=1000,
        r_EU=0.004,
        N_AS0=510,
        r_AS=0.0055,
        m_EU_AS=0.000096,
        T_AF=220000//25, 
        T_B=140000//25, 
        T_EU_AS=21200//25, 
        ops=[],
        infoFields=[],
        scale=1
        ):
        '''Counting **backward in time**, this model evolves a population for ``T0``
        generations (required parameter). The ancient population ``A`` started at
        size ``N_A`` and expanded at ``T_AF`` generations from now, to pop ``AF``
        with size ``N_AF``. Pop ``B`` split from pop ``AF`` at ``T_B`` generations
        from now, with size ``N_B``; Pop ``AF`` remains as ``N_AF`` individuals. 
        Pop ``EU`` and  ``AS`` split from pop ``B`` at ``T_EU_AS`` generations
        from now; with size ``N_EU0`` individuals and ``N_ASO`` individuals,
        respectively. Pop ``EU`` grew exponentially with rate ``r_EU``; Pop
        ``AS`` grew exponentially with rate ``r_AS``. The ``YRI``, ``CEU`` and
        ``CHB`` samples are drawn from ``AF``, ``EU`` and ``AS`` populations
        respectively. Additional operators could be added to ``ops``. Information
        fields required by these operators should be passed to ``infoFields``. If 
        a scaling factor ``scale`` is specified, all population sizes and
        generation numbers will be divided by a factor of ``scale``.

        This model merges all subpopulations if it is applied to a population with
        multiple subpopulation.
        '''
        #
        if T0 < T_AF:
            raise ValueError('Length of evolution T0=%d should be more than T_AF=%d' % (T0, T_AF))
        # for python 2.x and 3.x compatibility
        scale = float(scale)
        MultiStageModel.__init__(self, [
            InstantChangeModel(
                T=int((T0-T_EU_AS)/scale),
                N0=(int(N_A/scale), 'Ancestral'),
                # change population size twice, one at T_AF, one at T_B
                G=[int((T0-T_AF)/scale), int((T0-T_B)/scale)],
                NG=[
                    (int(N_AF/scale), 'AF'), 
                    # at T_B, split to population B from subpopulation 1
                    (int(N_B/scale), 'B')]),
            #
            ExponentialGrowthModel(
                T=int(T_EU_AS/scale),
                # split B into EU and AS at the beginning of this
                # exponential growth stage
                N0 = [(int(N_EU0/scale), 'EU'), (int(N_AS0/scale), 'AS')],
                r=[r_EU*scale, r_AS*scale],
                infoFields='migrate_to',
                ops=Migrator(rate=[
                    [0, m_EU_AS],
                    [m_EU_AS, 0]
                    ])
                ),
            ], ops=ops, infoFields=infoFields
        )

