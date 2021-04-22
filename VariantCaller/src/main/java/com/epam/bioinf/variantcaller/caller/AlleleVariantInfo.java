package com.epam.bioinf.variantcaller.caller;

import com.epam.bioinf.variantcaller.caller.variant.VariantInfo;
import htsjdk.variant.variantcontext.Allele;

import java.util.Arrays;
import java.util.Collection;
import java.util.stream.Stream;


public class AlleleVariantInfo {
    Allele[] alleles;
    VariantInfo[] variantInfos;
    int capacity; // overall mem allocation
    int size; // current amount of elements

    AlleleVariantInfo(final Allele allele, final VariantInfo variantInfo) {
        size = 1;
        capacity = 1;
        alleles = new Allele[1];
        alleles[0] = allele;
        variantInfos = new VariantInfo[1];
        variantInfos[0] = variantInfo;
    }

    void put(final Allele allele, final VariantInfo variantInfo) {
        if (alleleExists(allele)) {
            return; // allele already present -> not adding
        }
        Allele[] tempAlleles;
        VariantInfo[] tempVariantInfos;
        if (capacity == size) {
            capacity *= 2;
            tempAlleles = new Allele[capacity];
            tempVariantInfos = new VariantInfo[capacity];
            System.arraycopy(alleles, 0, tempAlleles, 0, size);
            System.arraycopy(variantInfos, 0, tempVariantInfos, 0, size);
            alleles = tempAlleles;
            variantInfos = tempVariantInfos;
        }
        alleles[++size] = allele;
        variantInfos[++size] = variantInfo;
    }

    public VariantInfo get(final Allele allele) {
        var index = -1;
        for(int i = 0; i < alleles.length; i++) {
            if (alleles[i].equals(allele)) {
                index = i;
                break;
            }
        }
        if (index == -1) {
            return null;
        }
        return variantInfos[index];
    }

    private boolean alleleExists(final Allele allele) {
        for (Allele value : alleles) {
            if (value.equals(allele)) {
                return true;
            }
        }
        return false;
    }

    public Collection<VariantInfo> values() {
        return Arrays.asList(variantInfos);
    }
}
