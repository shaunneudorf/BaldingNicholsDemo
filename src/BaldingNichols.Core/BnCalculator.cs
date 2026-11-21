namespace BaldingNichols.Core;

/// <summary>
/// One draw of the unknown contributor's allele using the Balding-Nichols
/// sampling formula:  [x·θ + (1−θ)·p] / [1 + (n−1)·θ]
/// x = copies of this allele already seen, n = total alleles already seen.
/// </summary>
public sealed record Draw(int Step, string Allele, int X, int N, double Frequency, double Theta)
{
    public double Numerator => X * Theta + (1 - Theta) * Frequency;
    public double Denominator => 1 + (N - 1) * Theta;
    public double Term => Numerator / Denominator;
}

/// <summary>One ordering of the unknown's two alleles (e.g. 12 then 13).</summary>
public sealed record OrderResult(string Label, IReadOnlyList<Draw> Draws)
{
    public double Probability => Draws.Aggregate(1.0, (acc, d) => acc * d.Term);
}

/// <summary>Result for one locus.</summary>
public sealed record LocusResult(
    double Theta,
    IReadOnlyList<string> Victim,
    IReadOnlyList<string> Suspect,
    IReadOnlyList<string> Conditioning,
    IReadOnlyList<OrderResult> Orders)
{
    public bool IsHomozygous => string.Equals(Suspect[0], Suspect[1], StringComparison.OrdinalIgnoreCase);

    /// <summary>Pr(unknown has the suspect's genotype | conditioning alleles).</summary>
    public double GenotypeProbability => Orders.Sum(o => o.Probability);

    /// <summary>LR = Pr(E|Hp) / Pr(E|Hd), with Pr(E|Hp) = 1.</summary>
    public double LikelihoodRatio => 1.0 / GenotypeProbability;
}

/// <summary>Product-rule (θ = 0) result for the suspect's genotype.</summary>
public sealed record ProductRuleResult(IReadOnlyList<string> Suspect, double P1, double P2)
{
    public bool IsHomozygous => string.Equals(Suspect[0], Suspect[1], StringComparison.OrdinalIgnoreCase);
    public double Probability => IsHomozygous ? P1 * P1 : 2 * P1 * P2;
    public double LikelihoodRatio => 1.0 / Probability;
}

public static class BnCalculator
{
    private static readonly StringComparer Cmp = StringComparer.OrdinalIgnoreCase;

    /// <param name="theta">Co-ancestry coefficient θ (FST), 0 ≤ θ &lt; 1.</param>
    /// <param name="victim">Assumed contributor's alleles: empty (none) or two alleles.</param>
    /// <param name="suspect">Suspect's two alleles (same allele twice = homozygote).</param>
    /// <param name="frequencies">Population frequency for each suspect allele.</param>
    public static LocusResult Compute(
        double theta,
        IEnumerable<string> victim,
        IEnumerable<string> suspect,
        IReadOnlyDictionary<string, double> frequencies)
    {
        if (double.IsNaN(theta) || theta < 0 || theta >= 1)
            throw new ArgumentException("θ must be at least 0 and less than 1.");

        var v = victim.Select(Norm).Where(a => a.Length > 0).ToList();
        var s = suspect.Select(Norm).ToList();

        if (s.Count != 2 || s.Any(a => a.Length == 0))
            throw new ArgumentException("Enter both suspect alleles.");
        if (v.Count != 0 && v.Count != 2)
            throw new ArgumentException("Enter both victim alleles, or leave both blank.");

        var freq = Canonical(frequencies);
        foreach (var a in s.Distinct(Cmp))
        {
            if (!freq.TryGetValue(a, out var p) || double.IsNaN(p) || p <= 0 || p > 1)
                throw new ArgumentException($"Enter a frequency between 0 and 1 for allele {a}.");
        }

        // The alleles already "seen": the victim's and the suspect's.
        var pool = v.Concat(s).ToList();

        var orders = new List<OrderResult> { BuildOrder("Order 1", pool, s[0], s[1], freq, theta) };
        if (!Cmp.Equals(s[0], s[1]))
            orders.Add(BuildOrder("Order 2", pool, s[1], s[0], freq, theta));

        return new LocusResult(theta, v, s, pool, orders);
    }

    /// <summary>Plain product rule (θ = 0): 2pq for a heterozygote, p² for a homozygote.</summary>
    public static ProductRuleResult ProductRule(IEnumerable<string> suspect, IReadOnlyDictionary<string, double> frequencies)
    {
        var s = suspect.Select(Norm).ToList();
        if (s.Count != 2 || s.Any(a => a.Length == 0))
            throw new ArgumentException("Enter both suspect alleles.");

        var freq = Canonical(frequencies);
        double Get(string a) =>
            freq.TryGetValue(a, out var p) && !double.IsNaN(p) && p > 0 && p <= 1
                ? p
                : throw new ArgumentException($"Enter a frequency between 0 and 1 for allele {a}.");

        return new ProductRuleResult(s, Get(s[0]), Get(s[1]));
    }

    public static double ProductRuleProbability(IEnumerable<string> suspect, IReadOnlyDictionary<string, double> frequencies) =>
        ProductRule(suspect, frequencies).Probability;

    private static OrderResult BuildOrder(
        string label, List<string> pool, string first, string second,
        Dictionary<string, double> freq, double theta)
    {
        var draws = new List<Draw>();
        var drawn = new List<string>();
        foreach (var allele in new[] { first, second })
        {
            int x = pool.Count(a => Cmp.Equals(a, allele)) + drawn.Count(a => Cmp.Equals(a, allele));
            int n = pool.Count + drawn.Count;
            draws.Add(new Draw(drawn.Count + 1, allele, x, n, freq[allele], theta));
            drawn.Add(allele);
        }
        return new OrderResult(label, draws);
    }

    private static string Norm(string? a) => AlleleName.Canon(a);

    private static Dictionary<string, double> Canonical(IReadOnlyDictionary<string, double> frequencies)
    {
        var d = new Dictionary<string, double>(Cmp);
        foreach (var kv in frequencies) d[Norm(kv.Key)] = kv.Value;
        return d;
    }
}
