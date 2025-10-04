"""
Alias execution CLI for qxub.

Provides commands for executing aliases with optional overrides and command arguments.
"""
# pylint: disable=import-outside-toplevel

import click
from rich.console import Console

from .config_manager import config_manager


console = Console()


@click.command(name='alias')
@click.argument('alias_name')
@click.argument('command_args', nargs=-1)
@click.option('--cmd', help='Override the alias command')
@click.option('--name', help='Override job name')
@click.option('--queue', help='Override queue')
@click.option('--project', help='Override project')
@click.option('--resources', multiple=True, help='Override resource requests')
@click.option('--joblog', help='Override joblog path')
@click.option('--out', help='Override stdout path')
@click.option('--err', help='Override stderr path')
@click.option('--pre', help='Override pre-command')
@click.option('--post', help='Override post-command')
# Conda-specific overrides
@click.option('--env', help='Override conda environment')
# Module-specific overrides
@click.option('--mod', help='Override single module to load')
@click.option('--mods', help='Override multiple modules to load (comma-separated)')
# Singularity-specific overrides
@click.option('--sif', help='Override singularity container')
@click.option('--bind', multiple=True, help='Override bind paths')
@click.option('--env-var', multiple=True, help='Override environment variables')
@click.pass_context
def alias_cli(ctx, alias_name: str, command_args: tuple, **overrides):  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    """Execute an alias with optional overrides and command arguments.

    Examples:
        qxub alias dvc_push
        qxub alias train_model --queue gpuvolta
        qxub alias analysis -- input.bam output.bam
        qxub alias quick_task --cmd "python script.py"
        qxub alias module_task --mod python3
        qxub alias multi_task --mods "python3,samtools,gcc"
    """
    # Get the alias definition
    alias_def = config_manager.get_alias(alias_name)
    if not alias_def:
        click.echo(f"❌ Alias '{alias_name}' not found")
        click.echo("Available aliases:")
        for alias in config_manager.list_aliases():
            click.echo(f"  • {alias}")
        raise click.Abort()

    # Check if alias has a subcommand
    subcommand = alias_def.get('subcommand')
    if not subcommand:
        click.echo(f"❌ Alias '{alias_name}' missing subcommand")
        raise click.Abort()

    if subcommand not in ['conda', 'module', 'sing']:
        click.echo(f"❌ Invalid subcommand in alias '{alias_name}': {subcommand}")
        raise click.Abort()

    # Prepare the command
    cmd = overrides.get('cmd') or alias_def.get('cmd')
    if not cmd and not command_args:
        click.echo(f"❌ No command specified for alias '{alias_name}'")
        click.echo("Either define 'cmd' in the alias or provide command arguments")
        raise click.Abort()

    # Append command arguments if provided
    if command_args:
        if cmd:
            cmd = f"{cmd} {' '.join(command_args)}"
        else:
            cmd = ' '.join(command_args)

    # Remove None values from overrides
    clean_overrides = {k: v for k, v in overrides.items() if v is not None}

    # Validate mutually exclusive module options
    if 'mod' in clean_overrides and 'mods' in clean_overrides:
        click.echo("❌ Cannot specify both --mod and --mods")
        raise click.Abort()

    # Handle module options
    if 'mod' in clean_overrides:
        # Single module - convert to list for internal consistency
        clean_overrides['mod'] = [clean_overrides['mod']]
    elif 'mods' in clean_overrides:
        # Multiple modules - split by comma and clean up
        modules = [m.strip() for m in clean_overrides['mods'].split(',') if m.strip()]
        clean_overrides['mod'] = modules
        del clean_overrides['mods']  # Remove mods key, use mod for consistency

    # Convert multiple values to lists
    if 'resources' in clean_overrides and clean_overrides['resources']:
        clean_overrides['resources'] = list(clean_overrides['resources'])
    if 'bind' in clean_overrides and clean_overrides['bind']:
        clean_overrides['bind'] = list(clean_overrides['bind'])
    if 'env_var' in clean_overrides and clean_overrides['env_var']:
        clean_overrides['env_var'] = list(clean_overrides['env_var'])

    # Get global options from parent context
    parent_ctx = ctx.parent
    global_options = {}
    if parent_ctx:
        global_options = {
            'dry': parent_ctx.params.get('dry', False),
            'quiet': parent_ctx.params.get('quiet', False),
            'verbose': parent_ctx.params.get('verbose', False),
        }

    # Resolve final options
    resolved_options = config_manager.resolve_options(clean_overrides, alias_name)

    # Extract subcommand-specific options
    subcommand_options = resolved_options.get(subcommand, {})

    # Handle subcommand-specific overrides
    if subcommand == 'conda' and 'env' in clean_overrides:
        subcommand_options['env'] = clean_overrides['env']
    elif subcommand == 'module' and 'mod' in clean_overrides and clean_overrides['mod']:
        # Only override if mod is not empty
        subcommand_options['mod'] = clean_overrides['mod']
    elif subcommand == 'sing':
        if 'sif' in clean_overrides:
            subcommand_options['sif'] = clean_overrides['sif']
        if 'bind' in clean_overrides and clean_overrides['bind']:
            # Only override if bind is not empty
            subcommand_options['bind'] = clean_overrides['bind']
        if 'env_var' in clean_overrides and clean_overrides['env_var']:
            # Only override if env_var is not empty
            subcommand_options['env'] = clean_overrides['env_var']

    # Import subcommands (avoid circular imports)
    from .conda import conda
    from .module import module
    from .sing import sing

    subcommand_map = {
        'conda': conda,
        'module': module,
        'sing': sing
    }

    # Prepare arguments for the subcommand
    subcommand_func = subcommand_map[subcommand]

    # Build argument list
    args = [cmd]  # The command is always the first argument

    # Add main options
    if resolved_options.get('name'):
        args.extend(['--name', resolved_options['name']])
    if resolved_options.get('queue'):
        args.extend(['--queue', resolved_options['queue']])
    if resolved_options.get('project'):
        args.extend(['--project', resolved_options['project']])
    if resolved_options.get('joblog'):
        args.extend(['--joblog', resolved_options['joblog']])
    if resolved_options.get('out'):
        args.extend(['--out', resolved_options['out']])
    if resolved_options.get('err'):
        args.extend(['--err', resolved_options['err']])
    if resolved_options.get('pre'):
        args.extend(['--pre', resolved_options['pre']])
    if resolved_options.get('post'):
        args.extend(['--post', resolved_options['post']])

    # Add resources
    resources = resolved_options.get('resources', [])
    if isinstance(resources, str):
        resources = [resources]
    for resource in resources:
        args.extend(['--resources', resource])

    # Add subcommand-specific options
    if subcommand == 'conda':
        if subcommand_options.get('env'):
            args.extend(['--env', subcommand_options['env']])
    elif subcommand == 'module':
        modules = subcommand_options.get('mod', [])
        if isinstance(modules, str):
            modules = [modules]
        for mod in modules:
            args.extend(['--mod', mod])
    elif subcommand == 'sing':
        if subcommand_options.get('sif'):
            args.extend(['--sif', subcommand_options['sif']])

        bind_paths = subcommand_options.get('bind', [])
        if isinstance(bind_paths, str):
            bind_paths = [bind_paths]
        for bind_path in bind_paths:
            args.extend(['--bind', bind_path])

        env_vars = subcommand_options.get('env', [])
        if isinstance(env_vars, str):
            env_vars = [env_vars]
        for env_var in env_vars:
            args.extend(['--env', env_var])

    # Show what we're about to execute if verbose or dry run
    if global_options.get('verbose') or global_options.get('dry'):
        click.echo(f"🚀 Executing alias: {alias_name}")
        click.echo(f"📋 Subcommand: qxub {subcommand}")
        click.echo(f"💻 Command: {cmd}")
        if global_options.get('dry'):
            click.echo(f"🔍 Full args: {args}")
            return

    # Create a new context for the subcommand
    try:
        # We need to invoke the subcommand with proper context
        with subcommand_func.make_context(subcommand, args, parent=ctx) as sub_ctx:
            # Apply global options to the new context
            sub_ctx.parent.params.update(global_options)

            # Invoke the subcommand
            subcommand_func.invoke(sub_ctx)

    except click.ClickException:
        # Re-raise click exceptions
        raise
    except Exception as e:
        click.echo(f"❌ Error executing alias '{alias_name}': {e}")
        raise click.Abort()


# Test command for aliases
@click.command(name='alias-test')
@click.argument('alias_name')
def alias_test_cli(alias_name: str):
    """Test an alias without executing it (dry run)."""
    # Get the alias definition
    alias_def = config_manager.get_alias(alias_name)
    if not alias_def:
        click.echo(f"❌ Alias '{alias_name}' not found")
        raise click.Abort()

    click.echo(f"🧪 Testing alias: {alias_name}")

    # Show resolved configuration
    resolved_options = config_manager.resolve_options({}, alias_name)

    click.echo("📋 Resolved configuration:")
    for key, value in resolved_options.items():
        if isinstance(value, dict):
            click.echo(f"  {key}:")
            for subkey, subvalue in value.items():
                click.echo(f"    {subkey}: {subvalue}")
        else:
            click.echo(f"  {key}: {value}")

    # Validate the alias
    subcommand = alias_def.get('subcommand')
    cmd = alias_def.get('cmd')

    click.echo("✅ Alias validation:")
    click.echo(f"  • Subcommand: {subcommand}")
    click.echo(f"  • Command: {cmd or '(none - requires command args)'}")

    if subcommand == 'sing':
        sif = resolved_options.get('sing', {}).get('sif')
        if not sif:
            click.echo("  ⚠️  Warning: No Singularity container specified")
        else:
            click.echo(f"  • Container: {sif}")

    click.echo("🎉 Alias test completed")
