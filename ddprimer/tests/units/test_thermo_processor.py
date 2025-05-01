#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the thermodynamic processor module.
Tests the flexible backend system that can switch between NUPACK and ViennaRNA.
"""

import pandas as pd
import pytest
from unittest import mock

# Import the module to test
from ...core import ThermoProcessor
from ...config import Config


class TestThermoProcessor:
    """Test suite for ThermoProcessor class."""

    def setup_method(self):
        """Set up test data before each test."""
        # Mock configuration settings
        self.original_thermo_backend = getattr(Config, 'THERMO_BACKEND', None)
        self.original_thermo_auto_fallback = getattr(Config, 'THERMO_AUTO_FALLBACK', None)
        self.original_thermo_temperature = getattr(Config, 'THERMO_TEMPERATURE', None)
        self.original_thermo_sodium = getattr(Config, 'THERMO_SODIUM', None)
        self.original_thermo_magnesium = getattr(Config, 'THERMO_MAGNESIUM', None)
        self.original_show_progress = getattr(Config, 'SHOW_PROGRESS', None)
        
        # Set test config values
        Config.THERMO_BACKEND = 'vienna'  # Default to ViennaRNA for tests
        Config.THERMO_AUTO_FALLBACK = True
        Config.THERMO_TEMPERATURE = 37
        Config.THERMO_SODIUM = 0.05
        Config.THERMO_MAGNESIUM = 0.002
        Config.SHOW_PROGRESS = False
        
        # Reset cached backend before each test
        ThermoProcessor._backend = None
        
        # Sample sequences for testing
        self.test_sequence = "ATCGATCGATCGTAGCTAGCTAGC"
        self.test_sequences = [
            "ATCGATCGATCG",
            "GCTAGCTAGCTA",
            "NNNNNNNNNNNN",  # All Ns
            ""  # Empty sequence
        ]
        
        # Sample pandas Series for testing
        self.test_series = pd.Series(self.test_sequences)

    def teardown_method(self):
        """Restore configuration after each test."""
        if self.original_thermo_backend is not None:
            Config.THERMO_BACKEND = self.original_thermo_backend
        if self.original_thermo_auto_fallback is not None:
            Config.THERMO_AUTO_FALLBACK = self.original_thermo_auto_fallback
        if self.original_thermo_temperature is not None:
            Config.THERMO_TEMPERATURE = self.original_thermo_temperature
        if self.original_thermo_sodium is not None:
            Config.THERMO_SODIUM = self.original_thermo_sodium
        if self.original_thermo_magnesium is not None:
            Config.THERMO_MAGNESIUM = self.original_thermo_magnesium
        if self.original_show_progress is not None:
            Config.SHOW_PROGRESS = self.original_show_progress
        
        # Reset cached backend after each test
        ThermoProcessor._backend = None

    @mock.patch('ddprimer.core.thermo_processor.ThermoProcessor.check_nupack_available')
    @mock.patch('ddprimer.core.thermo_processor.ThermoProcessor.check_vienna_available')
    def test_get_backend(self, mock_vienna_available, mock_nupack_available):
        """Test backend selection logic."""
        # Set up mocks
        mock_nupack_available.return_value = True
        mock_vienna_available.return_value = True
        
        # Test default to configured backend (vienna)
        Config.THERMO_BACKEND = 'vienna'
        backend = ThermoProcessor.get_backend()
        assert backend == 'vienna'
        
        # Test using NUPACK when configured
        Config.THERMO_BACKEND = 'nupack'
        ThermoProcessor._backend = None  # Reset cached value
        backend = ThermoProcessor.get_backend()
        assert backend == 'nupack'
        
        # Test fallback when configured backend is not available
        Config.THERMO_BACKEND = 'vienna'
        mock_vienna_available.return_value = False
        ThermoProcessor._backend = None  # Reset cached value
        backend = ThermoProcessor.get_backend()
        assert backend == 'nupack'  # Should fall back to NUPACK
        
        # Test no fallback when auto_fallback is disabled
        Config.THERMO_AUTO_FALLBACK = False
        ThermoProcessor._backend = None  # Reset cached value
        backend = ThermoProcessor.get_backend()
        assert backend is None  # Should not fall back
        
        # Test when no backend is available
        mock_nupack_available.return_value = False
        mock_vienna_available.return_value = False
        Config.THERMO_AUTO_FALLBACK = True
        ThermoProcessor._backend = None  # Reset cached value
        backend = ThermoProcessor.get_backend()
        assert backend is None

    def test_set_backend(self):
        """Test explicitly setting the backend."""
        # Mock the availability checks
        with mock.patch('ddprimer.core.thermo_processor.ThermoProcessor.check_nupack_available') as mock_nupack, \
             mock.patch('ddprimer.core.thermo_processor.ThermoProcessor.check_vienna_available') as mock_vienna:
            
            # Both backends available
            mock_nupack.return_value = True
            mock_vienna.return_value = True
            
            # Set to NUPACK
            result = ThermoProcessor.set_backend('nupack')
            assert result is True
            assert ThermoProcessor._backend == 'nupack'
            
            # Set to ViennaRNA
            result = ThermoProcessor.set_backend('vienna')
            assert result is True
            assert ThermoProcessor._backend == 'vienna'
            
            # Try to set to non-existent backend
            result = ThermoProcessor.set_backend('nonexistent')
            assert result is False
            assert ThermoProcessor._backend == 'vienna'  # Should not change
            
            # Try to set to unavailable backend
            mock_nupack.return_value = False
            result = ThermoProcessor.set_backend('nupack')
            assert result is False
            assert ThermoProcessor._backend == 'vienna'  # Should not change

    @mock.patch('ddprimer.core.nupack_processor.NupackProcessor.calc_deltaG')
    @mock.patch('ddprimer.core.vienna_processor.ViennaRNAProcessor.calc_deltaG')
    @mock.patch('ddprimer.core.thermo_processor.ThermoProcessor.get_backend')
    def test_calc_deltaG(self, mock_get_backend, mock_vienna_calc, mock_nupack_calc):
        """Test deltaG calculation with different backends."""
        # Set up mock return values
        mock_vienna_calc.return_value = -10.5
        mock_nupack_calc.return_value = -9.8
        
        # Test ViennaRNA backend
        mock_get_backend.return_value = 'vienna'
        result = ThermoProcessor.calc_deltaG(self.test_sequence)
        assert result == -10.5
        assert mock_vienna_calc.called
        assert not mock_nupack_calc.called
        
        # Reset mocks
        mock_vienna_calc.reset_mock()
        mock_nupack_calc.reset_mock()
        
        # Test NUPACK backend
        mock_get_backend.return_value = 'nupack'
        result = ThermoProcessor.calc_deltaG(self.test_sequence)
        assert result == -9.8
        assert not mock_vienna_calc.called
        assert mock_nupack_calc.called
        
        # Reset mocks
        mock_vienna_calc.reset_mock()
        mock_nupack_calc.reset_mock()
        
        # Test when no backend is available
        mock_get_backend.return_value = None
        result = ThermoProcessor.calc_deltaG(self.test_sequence)
        assert result is None
        assert not mock_vienna_calc.called
        assert not mock_nupack_calc.called

    @mock.patch('ddprimer.core.nupack_processor.NupackProcessor.calc_deltaG_batch')
    @mock.patch('ddprimer.core.vienna_processor.ViennaRNAProcessor.calc_deltaG_batch')
    @mock.patch('ddprimer.core.thermo_processor.ThermoProcessor.get_backend')
    def test_calc_deltaG_batch(self, mock_get_backend, mock_vienna_batch, mock_nupack_batch):
        """Test batch deltaG calculation with different backends."""
        # Set up mock return values
        mock_vienna_batch.return_value = [-10.5, -9.2, None, None]
        mock_nupack_batch.return_value = [-9.8, -8.5, None, None]
        
        # Test ViennaRNA backend
        mock_get_backend.return_value = 'vienna'
        result = ThermoProcessor.calc_deltaG_batch(self.test_sequences)
        assert result == [-10.5, -9.2, None, None]
        assert mock_vienna_batch.called
        assert not mock_nupack_batch.called
        
        # Reset mocks
        mock_vienna_batch.reset_mock()
        mock_nupack_batch.reset_mock()
        
        # Test NUPACK backend
        mock_get_backend.return_value = 'nupack'
        result = ThermoProcessor.calc_deltaG_batch(self.test_sequences)
        assert result == [-9.8, -8.5, None, None]
        assert not mock_vienna_batch.called
        assert mock_nupack_batch.called
        
        # Reset mocks
        mock_vienna_batch.reset_mock()
        mock_nupack_batch.reset_mock()
        
        # Test when no backend is available
        mock_get_backend.return_value = None
        result = ThermoProcessor.calc_deltaG_batch(self.test_sequences)
        assert result == [None, None, None, None]
        assert not mock_vienna_batch.called
        assert not mock_nupack_batch.called
        
        # Test with custom description
        mock_get_backend.return_value = 'vienna'
        result = ThermoProcessor.calc_deltaG_batch(self.test_sequences, "Custom description")
        assert mock_vienna_batch.call_args[0][1] == "Custom description"

    @mock.patch('ddprimer.core.nupack_processor.NupackProcessor.process_deltaG')
    @mock.patch('ddprimer.core.vienna_processor.ViennaRNAProcessor.process_deltaG')
    @mock.patch('ddprimer.core.thermo_processor.ThermoProcessor.get_backend')
    def test_process_deltaG(self, mock_get_backend, mock_vienna_process, mock_nupack_process):
        """Test pandas Series processing with different backends."""
        # Set up mock return values
        mock_vienna_process.return_value = pd.Series([-10.5, -9.2, None, None])
        mock_nupack_process.return_value = pd.Series([-9.8, -8.5, None, None])
        
        # Test ViennaRNA backend
        mock_get_backend.return_value = 'vienna'
        result = ThermoProcessor.process_deltaG(self.test_series)
        pd.testing.assert_series_equal(result, pd.Series([-10.5, -9.2, None, None]))
        assert mock_vienna_process.called
        assert not mock_nupack_process.called
        
        # Reset mocks
        mock_vienna_process.reset_mock()
        mock_nupack_process.reset_mock()
        
        # Test NUPACK backend
        mock_get_backend.return_value = 'nupack'
        result = ThermoProcessor.process_deltaG(self.test_series)
        pd.testing.assert_series_equal(result, pd.Series([-9.8, -8.5, None, None]))
        assert not mock_vienna_process.called
        assert mock_nupack_process.called
        
        # Reset mocks
        mock_vienna_process.reset_mock()
        mock_nupack_process.reset_mock()
        
        # Test when no backend is available
        mock_get_backend.return_value = None
        result = ThermoProcessor.process_deltaG(self.test_series)
        pd.testing.assert_series_equal(result, pd.Series([None, None, None, None], index=self.test_series.index))
        assert not mock_vienna_process.called
        assert not mock_nupack_process.called
        
        # Test with custom description
        mock_get_backend.return_value = 'vienna'
        result = ThermoProcessor.process_deltaG(self.test_series, "Custom description")
        assert mock_vienna_process.call_args[0][1] == "Custom description"