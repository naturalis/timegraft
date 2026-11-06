"""Tools that are not packaged in conda channels, built from source once."""

# pinned commit of https://github.com/tothuhien/lsd2 (reports v.2.4.4; the latest tag is older)
LSD2_COMMIT = config.get("lsd2_commit", "c61110f3a4fa05325b45c97b2134792ff9d55d4c")
LSD2_BIN = f"resources/bin/lsd2-{LSD2_COMMIT[:8]}"
USE_LSD2 = config["markers"][config["marker"]].get("dating", {}).get("method", "mpl") == "lsd2"


rule build_lsd2:
    output:
        LSD2_BIN,
    params:
        commit=LSD2_COMMIT,
    log:
        f"{OUT}/logs/build_lsd2.log",
    conda:
        "../envs/lsd2.yaml"
    shell:
        """
        (
        tmp=$(mktemp -d)
        git clone --quiet https://github.com/tothuhien/lsd2.git "$tmp/lsd2"
        git -C "$tmp/lsd2" checkout --quiet {params.commit}
        make -C "$tmp/lsd2/src"
        mkdir -p $(dirname {output})
        cp "$tmp/lsd2/src/lsd2" {output}
        rm -rf "$tmp"
        ) > {log} 2>&1
        """
