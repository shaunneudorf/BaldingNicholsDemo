using System.Globalization;

namespace BaldingNichols.Core;

/// <summary>Normalises allele designations so "09.3", "9.30" and "9.3" all match.</summary>
public static class AlleleName
{
    public static string Canon(string? allele)
    {
        var a = (allele ?? "").Trim();
        return double.TryParse(a, NumberStyles.Float, CultureInfo.InvariantCulture, out var d)
               && !double.IsNaN(d) && !double.IsInfinity(d)
            ? d.ToString("0.########", CultureInfo.InvariantCulture)
            : a.ToUpperInvariant();
    }
}

public enum FrequencyMethod { BaldingNichols, ProductRule }

/// <summary>
/// A frequency looked up from the database, with the numbers behind it so the page can typeset the working.
/// Count is x (copies of the allele in the database), K the number of observed designations, N the number of alleles.
/// </summary>
public sealed record FrequencyLookup(
    double Value, bool Observed, string Explanation,
    FrequencyMethod Method, int Count, int K, int N);

/// <summary>Allele frequencies for one locus, plus the database size N and observed-allele count k.</summary>
public sealed class LocusFrequencies
{
    private readonly Dictionary<string, double> _byAllele = new(StringComparer.Ordinal);

    public LocusFrequencies(string name, int alleleCount, IEnumerable<KeyValuePair<string, double>> frequencies)
    {
        Name = name;
        N = alleleCount;
        Observed = frequencies.Where(kv => kv.Value > 0).ToList();
        foreach (var kv in Observed)
            _byAllele[AlleleName.Canon(kv.Key)] = kv.Value;

        if (N <= 0) throw new FormatException($"Locus {name} needs an N greater than 0.");
        if (K == 0) throw new FormatException($"Locus {name} has no observed alleles.");
    }

    public string Name { get; }

    /// <summary>N: the number of alleles in the database for this locus.</summary>
    public int N { get; }

    /// <summary>Alleles with a non-zero frequency, in file order.</summary>
    public IReadOnlyList<KeyValuePair<string, double>> Observed { get; }

    /// <summary>k: the number of allele designations with non-zero observations.</summary>
    public int K => Observed.Count;

    /// <summary>x: copies of the allele in the database, recovered as frequency × N (0 if not observed).</summary>
    public int Count(string allele) =>
        _byAllele.TryGetValue(AlleleName.Canon(allele), out var p) ? Math.Max(1, (int)Math.Round(p * N)) : 0;

    /// <summary>Balding-Nichols frequency of an unseen allele: (1/k) / (N+1), the x = 0 case of the general formula.</summary>
    public double UnseenBaldingNicholsFrequency => (1.0 / K) / (N + 1);

    /// <summary>
    /// Product-rule frequency of an unseen allele: 3/2N, where 2N is the number of alleles in the database.
    /// (In this file N is already the allele count, so this is 3/N.) Change the 3 here to use another minimum.
    /// </summary>
    public double UnseenProductRuleFrequency => ProductRuleUnseenNumerator / N;

    public const double ProductRuleUnseenNumerator = 3;

    /// <summary>Balding-Nichols frequency: (x + 1/k) / (N + 1), for observed and unseen alleles alike.</summary>
    public FrequencyLookup LookupBaldingNichols(string allele)
    {
        var x = Count(allele);
        var value = (x + 1.0 / K) / (N + 1);
        var text = $"({x} + 1/{K}) ÷ ({N} + 1)";
        return new FrequencyLookup(value, x > 0, x > 0 ? text : $"Not observed: {text}",
            FrequencyMethod.BaldingNichols, x, K, N);
    }

    /// <summary>Product-rule frequency: the database frequency, or 3/2N if the allele was not observed.</summary>
    public FrequencyLookup LookupProductRule(string allele)
    {
        if (_byAllele.TryGetValue(AlleleName.Canon(allele), out var p))
            return new FrequencyLookup(p, true, "Database frequency.", FrequencyMethod.ProductRule, Count(allele), K, N);

        return new FrequencyLookup(UnseenProductRuleFrequency, false, $"Not observed: 3 ÷ 2N = 3 ÷ {N}.",
            FrequencyMethod.ProductRule, 0, K, N);
    }
}

/// <summary>Allele frequency tables read from a CSV file.</summary>
public sealed class AlleleDatabase
{
    private AlleleDatabase(IReadOnlyList<LocusFrequencies> loci) => Loci = loci;

    public IReadOnlyList<LocusFrequencies> Loci { get; }

    public LocusFrequencies? Find(string? locus) =>
        Loci.FirstOrDefault(l => string.Equals(l.Name, locus?.Trim(), StringComparison.OrdinalIgnoreCase));

    /// <summary>
    /// Reads a CSV whose first column is the allele designation and whose other columns are loci.
    /// Each cell is a frequency. A row labelled "N" gives the number of alleles in the database for each locus.
    /// </summary>
    public static AlleleDatabase ParseCsv(string csv)
    {
        var lines = csv.TrimStart('\uFEFF')
            .Split('\n')
            .Select(l => l.TrimEnd('\r'))
            .Where(l => l.Trim().Length > 0)
            .ToList();

        if (lines.Count < 2) throw new FormatException("The allele database has no data rows.");

        var names = lines[0].Split(',').Skip(1).Select(c => c.Trim()).ToArray();
        var freqs = names.Select(_ => new List<KeyValuePair<string, double>>()).ToArray();
        var counts = new int[names.Length];
        var sawN = false;

        foreach (var line in lines.Skip(1))
        {
            var cells = line.Split(',');
            var label = cells[0].Trim();
            var isN = label.Equals("N", StringComparison.OrdinalIgnoreCase);
            sawN |= isN;

            for (int j = 0; j < names.Length; j++)
            {
                var text = j + 1 < cells.Length ? cells[j + 1].Trim() : "";
                var value = text.Length == 0 ? 0 : ParseNumber(text, label, names[j]);
                if (isN) counts[j] = (int)Math.Round(value);
                else freqs[j].Add(new KeyValuePair<string, double>(label, value));
            }
        }

        if (!sawN) throw new FormatException("The allele database needs a row labelled N with the number of alleles for each locus.");

        var loci = names.Select((name, j) => new LocusFrequencies(name, counts[j], freqs[j])).ToList();
        return new AlleleDatabase(loci);
    }

    private static double ParseNumber(string text, string row, string locus) =>
        double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out var d)
            ? d
            : throw new FormatException($"Could not read \"{text}\" for allele {row} at {locus}.");
}
