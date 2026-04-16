import argparse

def generate_distance_commands(chain_pair):
    """
    Generate ChimeraX distance commands for a specified chain pair.
    """
    residues = [16, 80, 224]
    commands = []
    
    if len(chain_pair) != 2:
        raise ValueError("Each chain pair must have exactly two chain identifiers.")
    
    chain1, chain2 = chain_pair
    for residue in residues:
        commands.append(f"distance #8/{chain1}:{residue}@CA #8/{chain2}:{residue}@CA")
    
    return commands

def main():
    parser = argparse.ArgumentParser(description="Generate ChimeraX commands for measuring distances between residues in specified chain pairs.")
    parser.add_argument("-c1", nargs=2, required=True, metavar=("CHAIN1", "CHAIN2"),
                        help="Specify the first chain pair (e.g., -c1 CQ CL).")
    parser.add_argument("-c2", nargs=2, required=False, metavar=("CHAIN1", "CHAIN2"),
                        help="Specify the second chain pair (e.g., -c2 CL CK).")
    
    args = parser.parse_args()
    
    try:
        # Generate commands for -c1
        commands = generate_distance_commands(args.c1)
        print("\n".join(commands))
        
        # Generate commands for -c2 if provided
        if args.c2:
            commands = generate_distance_commands(args.c2)
            print("\n".join(commands))
    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
