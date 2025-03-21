#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UI utilities for ddPrimer pipeline.
"""

import sys
import time
import threading


class UIUtils:
    """UI-related utility functions."""
    
    @staticmethod
    def spinner_task(stop_event):
        """
        Print a rotating spinner while stop_event is not set.
        Ensures proper cleanup when stopping.
        
        Args:
            stop_event (threading.Event): Event to signal when to stop the spinner
        """
        spinner = ['|', '/', '-', '\\']
        idx = 0
        sys.stdout.write(spinner[idx])
        sys.stdout.flush()
        while not stop_event.is_set():
            time.sleep(0.1)
            idx = (idx + 1) % 4
            sys.stdout.write("\b" + spinner[idx])
            sys.stdout.flush()
        
        # Make sure we clear the spinner character
        sys.stdout.write("\b \b")  # backspace, space, backspace to fully erase
        sys.stdout.flush()
    
    @staticmethod
    def start_spinner():
        """
        Start a spinner in a separate thread.
        
        Returns:
            tuple: (stop_event, spinner_thread) to stop the spinner later
        """
        stop_event = threading.Event()
        spinner_thread = threading.Thread(target=UIUtils.spinner_task, args=(stop_event,))
        spinner_thread.daemon = True
        spinner_thread.start()
        return stop_event, spinner_thread
    
    @staticmethod
    def stop_spinner(stop_event, spinner_thread):
        """
        Stop a spinner that was started with start_spinner.
        
        Args:
            stop_event (threading.Event): Event to signal when to stop the spinner
            spinner_thread (threading.Thread): Thread running the spinner
        """
        stop_event.set()
        spinner_thread.join()
