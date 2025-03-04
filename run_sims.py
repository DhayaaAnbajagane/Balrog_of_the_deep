#!/usr/bin/env python

#copied directly from https://github.com/beckermr/misc/blob/68b42ea226913b97236ce5145378c6cf6fe05b8b/bin/run-simple-des-y3-sim
import logging
import sys

import click
import yaml


@click.group()
def cli():
    """Run simple DES Y3 end-to-emd simulations."""
    pass


@cli.command()
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to simulate as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--config-file', type=str, required=True,
              help='the YAML config file')
def initialize(tilename, bands, output_desdata, config_file):
    
    from initializing import initialize_files

    with open(config_file, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.Loader)
    initialize_files(
        tilename=tilename,
        bands=[b for b in bands],
        output_desdata=output_desdata,
        config=config)
    
@cli.command()
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to prep for as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--n-files', type=int, default=None,
              help='number of SE images to keep - useful for testing')
@click.option('--config-file', type=str, required=True,
              help='the YAML config file')
def prep(tilename, bands, output_desdata, n_files, config_file):
    """Prepare a tile for simulating."""
    
    from band_infoing import make_band_info
    
    with open(config_file, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.Loader)
    make_band_info(
        tilename=tilename,
        bands=[b for b in bands],
        output_meds_dir=output_desdata,
        config=config['survey_kws'],
        n_files=n_files)

@cli.command()
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to simulate as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--seed', type=int, required=True,
              help='the base RNG seed')
@click.option('--config-file', type=str, required=True,
              help='the YAML config file')
def galsim(tilename, bands, output_desdata, seed, config_file):
    
    from simulating import End2EndSimulation
    
    
    with open(config_file, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.Loader)
    sim = End2EndSimulation(
        seed=seed,
        output_meds_dir=output_desdata,
        tilename=tilename,
        bands=[b for b in bands],
        gal_kws=config['gal_kws'],
        psf_kws=config['psf_kws'])
    sim.run()



@cli.command('swarp')
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to prep for as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--config-file', type=str, required=True,
              help='the YAML config file')
def swarp(tilename, bands, output_desdata, config_file):
    
    from coadding import MakeSwarpCoadds
    
    with open(config_file, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.Loader)
        
    coadd = MakeSwarpCoadds(
        output_meds_dir=output_desdata,
        tilename=tilename,
        bands=[b for b in bands],
        config=config)
    coadd.run()


@cli.command('source-extractor')
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to prep for as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--config-file', type=str, required=True,
              help='the YAML config file')
def source_extractor(tilename, bands, output_desdata, config_file):
    
    from srcextracting import MakeSrcExtractorCat
    
    with open(config_file, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.Loader)
        
    SrcExtractor = MakeSrcExtractorCat(
                                output_meds_dir=output_desdata,
                                tilename=tilename,
                                bands=[b for b in bands],
                                config=config)
    SrcExtractor.run()
    
    
@cli.command('subselect')
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to prep for as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--radius', type=float, required=True,
              help='the matching radius in arcmin')
def subselect(radius, tilename, bands, output_desdata):
    
    from subselecting import keep_injected_only
    
    keep_injected_only(radius   = radius, 
                       tilename = tilename, 
                       bands    = bands,
                       output_desdata = output_desdata, 
                       )


@cli.command()
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to make MEDS files for as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--config-file', type=str, required=True,
              help='the YAML config file')
@click.option('--meds-config-file', type=str, required=True,
              help='the YAML config file for MEDS making')
def meds(tilename, bands, output_desdata, config_file, meds_config_file):
    
    from medsing import make_meds_files
    
    with open(config_file, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.Loader)
    with open(meds_config_file, 'r') as fp:
        meds_config = yaml.load(fp, Loader=yaml.Loader)
    make_meds_files(
        tilename=tilename,
        bands=[b for b in bands],
        output_meds_dir=output_desdata,
        psf_kws=config['psf_kws'],
        meds_config=meds_config)


@cli.command()
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to simulate as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--seed', type=int, required=True,
              help='the base RNG seed')
@click.option('--metacal-config-file', type=str, required=True,
              help='the YAML config file for metacal run')
def metacal(tilename, bands, output_desdata, seed, metacal_config_file):
    
    from run_metacal import run_metacal
    
    
    with open(metacal_config_file, 'r') as fp:
        mcal_config = yaml.load(fp, Loader=yaml.Loader)
    run_metacal(
        tilename=tilename,
        output_meds_dir=output_desdata,
        bands=[b for b in bands],
        seed=seed,
        mcal_config=mcal_config)
    
    
@cli.command()
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to simulate as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--seed', type=int, required=True,
              help='the base RNG seed')
@click.option('--fitvd-config-file', type=str, required=True,
              help='the YAML config file for fitvd run')
@click.option('--shredx-config-file', type=str, required=True,
              help='the YAML config file for shredx run')
@click.option('--config-file', type=str, required=True,
              help='the YAML config file')
def fitvd(tilename, bands, output_desdata, seed, config_file, fitvd_config_file, shredx_config_file):
    
    from fitvding import MakeShredxCats, MakeFitvdCats
    
    
    with open(config_file, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.Loader)

#     Shredx = MakeShredxCats(output_meds_dir=output_desdata,
#                             tilename=tilename,
#                             bands=[b for b in bands],
#                             config=config,
#                             shredx_config_path = shredx_config_file)
#     Shredx.run()


    Fitvd = MakeFitvdCats(output_meds_dir=output_desdata,
                          tilename=tilename,
                          bands=[b for b in bands],
                          config=config,
                          fitvd_config_path = fitvd_config_file)

    Fitvd.run()


@cli.command()
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to simulate as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--config-file', type=str, required=True,
              help='the YAML config file')
def match(tilename, bands, output_desdata, config_file):
    
    from matching import match_catalogs
    
    with open(config_file, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.Loader)
    match_catalogs(
        tilename=tilename,
        bands=[b for b in bands],
        output_desdata=output_desdata,
        config=config)


@cli.command()
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to simulate as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--config-file', type=str, required=True,
              help='the YAML config file')
def finalize(tilename, bands, output_desdata, config_file):
    
    from finalizing import finalize_files
    
    with open(config_file, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.Loader)
    finalize_files(
        tilename=tilename,
        bands=[b for b in bands],
        output_desdata=output_desdata,
        config=config)

if __name__ == '__main__':
    cli()
