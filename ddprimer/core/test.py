def test_kmer_masking_standalone():
    """
    Standalone test to check k-mer masking impact.
    Run this to test if k-mer masking is the problem.
    """
    processor = Primer3Processor()
    
    # Sample template (replace with one from your actual data)
    sample_template = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    
    # Run debug
    result = processor.debug_kmer_masking_impact(sample_template, "test_sequence")
    
    print("\n=== K-MER MASKING TEST RESULTS ===")
    print(f"Masking enabled: {result.get('masking_enabled', False)}")
    print(f"K-mer files found: {result.get('kmer_files_found', 0)}")
    print(f"Template length: {result.get('template_length', 0)} bp")
    print(f"Primers without masking: {result.get('pairs_without_masking', 'N/A')}")
    print(f"Primers with masking: {result.get('pairs_with_masking', 'N/A')}")
    
    if result.get('pairs_without_masking', 0) > 0 and result.get('pairs_with_masking', 0) == 0:
        print("\n❌ K-MER MASKING IS BLOCKING ALL PRIMER DESIGN!")
        print("Recommendations:")
        print("1. Disable k-mer masking: export DDPRIMER_DISABLE_KMER_MASKING=true")
        print("2. Use different k-mer files for your organism")
        print("3. Check if your sequences match the k-mer organism")
    elif result.get('pairs_with_masking', 0) > 0:
        print("\n✅ K-mer masking is working correctly")
    else:
        print("\n⚠️  Neither version produced primers - check primer design parameters")