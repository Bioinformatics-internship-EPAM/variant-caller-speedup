package com.epam.bioinf.variantcaller.caller;

import com.epam.bioinf.variantcaller.caller.variant.VariantInfo;
import htsjdk.variant.variantcontext.Allele;

public class AlleleVariantInfo {
    Allele allele;
    VariantInfo variantInfo;

    AlleleVariantInfo(final Allele allele, final VariantInfo variantInfo) {
        this.allele = allele;
        this.variantInfo = variantInfo;
    }

    public VariantInfo get() {
        return variantInfo;
    }
}
